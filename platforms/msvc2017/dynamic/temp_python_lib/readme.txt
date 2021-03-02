This directory is a temporary home for the Python import library (.lib)
during the build. The MSVC project copies it here and deletes it when it's
done.

Why? We only want to maintain & distribute one MSVC project but it must
work with Python 2.5, 2.6 and 2.7. Since the paths to the include & library
files vary depending on the Python version, we have to update the paths 
programmatically via pre-build events within the MSVC project. 

Although it's not too hard to pass programmatically-generated options to the
compile step of a Visual Studio build, it's AFAICT impossible to do so with
the linker (see below). So instead we opted to hardcode a predictable path in
the MSVC build and simply ensure that the .lib file is there when the build
runs. This is a case of moving the mountain to Mohammed.

Since the .lib file is an import library, it's only used during the build and
not at runtime. Therefore we can delete it from its temporary location when 
the build is done.

This is a somewhat clumsy solution, but it's the best we could come up with.
Here's three other things that almost worked --

The @ option (to pass a compiler response file) is supported for link.exe.
It would be trivial to use code similar to write_python_include_path.py to 
write a file like "python_lib_path.rsp" that contains 
/LIBPATH:c:\python25\libs and then pass @python_lib_path.rsp to the linker.
However, the doc clearly states that the @ option for link.exe doesn't work
within MSVC.
ref: http://msdn.microsoft.com/en-us/library/4xdcbak7.aspx

The cl.exe /link option can pass options to link.exe, and the @ option works
for cl.exe in MSVC. So theoretically we could add to the
python_include_path.rsp file that's passed to cl.exe something like this:
   /link /LIBPATH:c:\python25\libs
   
But MSVC executes cl.exe with the /c option (== compile only). link.exe is
executed explicitly by MSVC, so cl.exe never gets the opportunity to pass 
the contents of the /link option to link.exe.

Finally, one can alter the build config on the fly via macros, but macros 
are not supported in MSVC 2008 Express Edition which we and (we expect) 
many of our users use.
