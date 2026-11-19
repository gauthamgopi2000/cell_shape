function [] = lammps_write(xyz, L, Nparticle, timestep,fname)
%%This function writes vertex positions in the LAMMPS dump format, so that
%%OVITO can read it

%XYZ_c = reshape(xyz, 3, [])';
XYZ_c = xyz(:,1:15);
natoms = Nparticle;


XYZ_c = XYZ_c';
if timestep-1 == 0
    fileID = fopen(fname,'w');
    else
    fileID = fopen(fname,'a');
end
fprintf(fileID,'ITEM: TIMESTEP\n');
fprintf(fileID, '%d\n', timestep);
fprintf(fileID, 'ITEM: NUMBER OF ATOMS\n');
fprintf(fileID, '%d\n', natoms);
fprintf(fileID, 'ITEM: BOX BOUNDS\n');
fprintf(fileID, '%d %d\n', [-L/2., L/2.]);
fprintf(fileID, '%d %d\n', [-L/2., L/2.]);
fprintf(fileID, '%d %d\n', [-L/2., L/2,]);
fprintf(fileID, 'ITEM: ATOMS id type x y z radius Nvertex A A0 L L0 mean_l mean_l0 phi mean_f_tangent\n');
fprintf(fileID,'%d %d %f %f %f %f %f %f %f %f %f %f %f %f %f\n', XYZ_c);


