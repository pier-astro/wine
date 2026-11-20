#include <Python.h>
#include "uniform.h"
#include "redwine_geometry.h"
#include <xsTypes.h>
#include <XSFunctions/Utilities/FunctionUtility.h>
#include <cstdlib>
#include <dlfcn.h>
#include <memory>
#include <string>

namespace {

struct PythonDeleter {
    void operator()(PyObject* object) const { Py_XDECREF(object); }
};

using PythonObject = std::unique_ptr<PyObject, PythonDeleter>;

std::string libraryDirectory() {
    Dl_info info{};
    if (!dladdr(reinterpret_cast<const void*>(&libraryDirectory), &info) || !info.dli_fname) {
        return {};
    }
    const std::unique_ptr<char, decltype(&std::free)> path(realpath(info.dli_fname, nullptr), &std::free);
    if (!path) {
        return {};
    }
    const std::string filename(path.get());
    return filename.substr(0, filename.rfind('/'));
}

// Resolve this at lmod time, while relative library paths are still valid.
const std::string nativeDirectory = libraryDirectory();

PyObject* importBackend() {
    if (!nativeDirectory.empty()) {
        // build/native sits two levels below the windkin Python package root.
        const std::string windkinRoot = nativeDirectory + "/../..";
        PyObject* sysPath = PySys_GetObject("path");  // Borrowed reference.
        if (!sysPath) {
            PyErr_SetString(PyExc_RuntimeError, "Python sys.path is unavailable");
            return nullptr;
        }
        for (const auto& directory : {windkinRoot + "/../../src", windkinRoot}) {
            const PythonObject entry(PyUnicode_DecodeFSDefault(directory.c_str()));
            if (!entry) {
                return nullptr;
            }
            const int present = PySequence_Contains(sysPath, entry.get());
            if (present < 0 || (present == 0 && PyList_Insert(sysPath, 0, entry.get()) < 0)) {
                return nullptr;
            }
        }
    }
    return PyImport_ImportModule("windkin.xspec");
}

class PythonLock {
public:
    PythonLock() {
        // XSPEC loads native models serially. Preserve a host interpreter
        // when the library is loaded from PyXspec.
        if (!Py_IsInitialized()) {
            Py_InitializeEx(0);
            PyEval_SaveThread();
        }
        state = PyGILState_Ensure();
    }

    ~PythonLock() { PyGILState_Release(state); }

private:
    PyGILState_STATE state;
};

void fail(const char* name) {
    if (PyErr_Occurred()) {
        PyErr_Print();
    }
    throw FunctionUtility::FunctionException(std::string("windkin: ") + name +
        " failed. Check the RED-WINE checkout and NumPy/SciPy for the build's Python.");
}

void evaluate(const char* name, bool additive, const RealArray& energy,
              const RealArray& params, RealArray& flux, RealArray& fluxErr) {
    const PythonLock lock;
    static PyObject* backend = nullptr;
    if (!backend) {
        // Keep the module alive for the lifetime of the embedded interpreter.
        backend = importBackend();
        if (!backend) fail(name);
    }

    const PythonObject function(PyObject_GetAttrString(backend, name));
    const size_t bins = energy.size() - 1;
    const PythonObject edges(PyTuple_New(energy.size()));
    const PythonObject parameters(PyTuple_New(params.size()));
    const PythonObject values(PyList_New(bins));
    if (!function || !edges || !parameters || !values) fail(name);
    for (size_t i = 0; i < energy.size(); ++i) {
        if (PyTuple_SetItem(edges.get(), i, PyFloat_FromDouble(energy[i])) < 0) fail(name);
    }
    for (size_t i = 0; i < params.size(); ++i) {
        if (PyTuple_SetItem(parameters.get(), i, PyFloat_FromDouble(params[i])) < 0) fail(name);
    }
    for (size_t i = 0; i < bins; ++i) {
        if (PyList_SetItem(values.get(), i, PyFloat_FromDouble(additive ? 0. : flux[i])) < 0) fail(name);
    }
    const PythonObject result(PyObject_CallFunctionObjArgs(
        function.get(), edges.get(), parameters.get(), values.get(), nullptr));
    if (!result || PyList_Size(values.get()) != static_cast<Py_ssize_t>(bins)) fail(name);
    flux.resize(bins);
    fluxErr.resize(0);
    for (size_t i = 0; i < bins; ++i) {
        flux[i] = PyFloat_AsDouble(PyList_GetItem(values.get(), i));
    }
    if (PyErr_Occurred()) fail(name);
}

}  // namespace

extern "C" void wuline(const RealArray& energy, const RealArray& params,
                        int, RealArray& flux, RealArray& fluxErr, const string&) {
    evaluate("wuline", true, energy, params, flux, fluxErr);
}

extern "C" void wgline(const RealArray& energy, const RealArray& params,
                        int, RealArray& flux, RealArray& fluxErr, const string&) {
    evaluate("wgline", true, energy, params, flux, fluxErr);
}

extern "C" void zwgline(const RealArray& energy, const RealArray& params,
                         int, RealArray& flux, RealArray& fluxErr, const string&) {
    evaluate("zwgline", true, energy, params, flux, fluxErr);
}

namespace {
void evaluateUniform(const RealArray& energy, const RealArray& params,
                     RealArray& flux, RealArray& fluxErr, bool los) {
    try {
        const std::vector<double> edges(std::begin(energy), std::end(energy));
        const std::vector<double> input(std::begin(flux), std::end(flux));
        std::array<double,6> parameters{{params[0],params[1],params[2],params[3],params[4],params[5]}};
        if (los) {
            if (!std::isfinite(params[2]) || !std::isfinite(params[3]) || !std::isfinite(params[4]))
                throw std::invalid_argument("wuconvlos_Cf: geometry must be finite");
            // XSPEC derivatives may probe beyond hard bounds. Continue at
            // the nearest boundary, as in RED-WINE's LOS model adapter.
            const double cover=windkin_uniform::clip(params[2],0.,1.);
            const double skew=windkin_uniform::clip(params[3],0.,1.);
            const double inclination=windkin_uniform::clip(params[4],0.,90.);
            const double radians=std::acos(-1.)/180;
            const auto geometry=redwine::losGeometry(cover,skew,inclination*radians);
            parameters[4]=inclination;
            parameters[2]=geometry.tin/radians;
            parameters[3]=geometry.tout/radians;
        }
        const auto result=windkin_uniform::convolve(edges,input,parameters);
        for (size_t i=0; i<result.size(); ++i) flux[i]=result[i];
        fluxErr.resize(0);
    } catch (const std::invalid_argument& error) {
        // Match the Python adapter's zero-flux behavior outside physical bounds.
        FunctionUtility::xsWrite(std::string("windkin: ")+error.what(),10);
        flux=0.;
        fluxErr.resize(0);
    } catch (const std::exception& error) {
        throw FunctionUtility::FunctionException(std::string("windkin: ")+error.what());
    }
}
} // namespace

extern "C" void wuconv(const RealArray& energy, const RealArray& params,
                        int, RealArray& flux, RealArray& fluxErr, const string&) {
    evaluateUniform(energy, params, flux, fluxErr, false);
}

extern "C" void wuconvlosCf(const RealArray& energy, const RealArray& params,
                            int, RealArray& flux, RealArray& fluxErr, const string&) {
    evaluateUniform(energy, params, flux, fluxErr, true);
}

extern "C" void wgconv(const RealArray& energy, const RealArray& params,
                        int, RealArray& flux, RealArray& fluxErr, const string&) {
    evaluate("wgconv", false, energy, params, flux, fluxErr);
}

extern "C" void wglinelosCf(const RealArray& energy, const RealArray& params,
                        int, RealArray& flux, RealArray& fluxErr, const string&) {
    evaluate("wglinelos_Cf", true, energy, params, flux, fluxErr);
}

extern "C" void zwglinelosCf(const RealArray& energy, const RealArray& params,
                         int, RealArray& flux, RealArray& fluxErr, const string&) {
    evaluate("zwglinelos_Cf", true, energy, params, flux, fluxErr);
}

extern "C" void wgconvlosCf(const RealArray& energy, const RealArray& params,
                        int, RealArray& flux, RealArray& fluxErr, const string&) {
    evaluate("wgconvlos_Cf", false, energy, params, flux, fluxErr);
}

extern "C" void wsect(const RealArray& energy, const RealArray& params,
                       int, RealArray& flux, RealArray& fluxErr, const string&) {
    evaluate("wsect", false, energy, params, flux, fluxErr);
}

extern "C" void wsectline(const RealArray& energy, const RealArray& params,
                           int, RealArray& flux, RealArray& fluxErr, const string&) {
    evaluate("wsectline", true, energy, params, flux, fluxErr);
}

extern "C" void zwsectline(const RealArray& energy, const RealArray& params,
                            int, RealArray& flux, RealArray& fluxErr, const string&) {
    evaluate("zwsectline", true, energy, params, flux, fluxErr);
}
