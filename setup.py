from setuptools import setup, Extension

module = Extension(
    "symnmf_ext",  # Module name to be imported in Python
    sources=["symnmfmodule.c", "symnmf.c"],  # C source files
    extra_compile_args=["-Wall", "-Wextra", "-Werror", "-pedantic-errors"]  # Compiler flags
)

setup(
    name="symnmf_ext",
    version="1.0",
    description="Python wrapper for C implementation of SymNMF",
    ext_modules=[module],
)
