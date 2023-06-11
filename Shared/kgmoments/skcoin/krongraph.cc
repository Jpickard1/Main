/**
 * @file ckrongraph.cc
 * A python module to generate kronecker graphs via coin-flipping.
 * @author David F. Gleich
 */

#include <Python.h>
#include <math.h>

#include "krongraph_coin_flip.cpp"


static PyObject*
ckrongraph_edges(PyObject *self, PyObject *args)
{
    int ok;
    PyObject* mat; // the initiator matrix
    int r;
    
    ok = PyArg_ParseTuple(args, "Oi:edges", &mat, &r);
    if (!ok) { return NULL; }
    
    if (r < 0) {
        PyErr_SetString(PyExc_ValueError, "Recursion value is negative");
        return NULL; 
    }
       
    mat = PySequence_Fast(mat, "argument must be iterable");
    if (!mat) { return NULL; }
    
    // check the size
    Py_ssize_t matlen = PySequence_Fast_GET_SIZE(mat);
    // check that matlen is a square
    Py_ssize_t sqrt_matlen = (Py_ssize_t)sqrt((double)matlen);
    if (sqrt_matlen*sqrt_matlen != matlen) {
        PyErr_SetString(PyExc_ValueError, "matrix argument must be square");
        Py_DECREF(mat);
        return NULL;
    }
    
    // handle a simple return
    if (r == 0) {
        Py_DECREF(mat);
        return Py_BuildValue("[(i,i)]",0,0);
    }
    
    // get matrix data
    std::vector<double> matdata(matlen);
    for (Py_ssize_t i=0; i<matlen; ++i) {
        PyObject *item = PySequence_Fast_GET_ITEM(mat, i);
        if (!PyFloat_Check(item)) {
            PyErr_SetString(PyExc_TypeError, "matrix item is not a float.");
            Py_DECREF(mat);
            return NULL;
        }
        matdata[i] = PyFloat_AsDouble(item);
    }
    
    Py_DECREF(mat);
    
    // setup the krongraph structure
    KronGraphCoinFlip kg(&matdata[0], (size_t)sqrt_matlen, (size_t)r);
    kg.generate_edges();
    
    PyObject* rval = PyList_New((Py_ssize_t)kg.edges.size());
    if (!rval) {
        return PyErr_NoMemory( );
    }
    for (size_t i=0; i<kg.edges.size(); ++i) {
        PyList_SET_ITEM(rval, (Py_ssize_t)i, 
            Py_BuildValue("(i,i)",kg.edges[i].first,kg.edges[i].second));
    }
    return rval;
}


PyMethodDef methods[] = {
    {"edges", ckrongraph_edges, METH_VARARGS, "Return a list of edges in a kronecker graph"},
    {NULL, NULL, 0, NULL}
};


PyMODINIT_FUNC 
initkrongraph()
{
    (void) Py_InitModule("krongraph", methods);   
}
