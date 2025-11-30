%MATLAB CODE %load epanet2.dll and .h
loadlibrary('epanet2', 'epanet2')
%display library's functions
libfunctions epanet2 ?full
%open network file, set output file eg. net1.inp
calllib('epanet2', 'ENopen', 'net1.inp','net1.txt','');
%run hydraulic simulation and save in output file
calllib('epanet2','ENsolveH');
%run quality simulations and save in output file (after ENsolveH)
calllib('epanet2','ENsolveQ');
%for methods with pointers (returning values in arguments) %setup a pointer with th
nnodes = libpointer('int32Ptr',0);
%calls the method
calllib('epanet2','ENgetcount', 0 , nnodes);
%returns the values stored in the pointer
get(nnodes,'Value');
%write a formated text report in output file
calllib('epanet2','ENreport');
%close toolkit
calllib('epanet2', 'ENclose');
%unload library
unloadlibrary epanet2