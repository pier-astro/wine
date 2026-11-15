#include <Python.h>
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

extern "C" void wuconv(const RealArray& energy, const RealArray& params,
                        int, RealArray& flux, RealArray& fluxErr, const string&) {
    evaluate("wuconv", false, energy, params, flux, fluxErr);
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
