function [clift,cdrag,cmom,fstat,clinv,clfs] = interpolation(aoa,cl,cd,cm,angle_of_attack,thickness,thick_prof,f_stat,cl_inv,cl_fs)
 % f_stat,cl_inv,cl_fs   
    clthick = zeros(6,1);
    cdthick = zeros(6,1);
    cmthick = zeros(6,1);
    fstat_thick = zeros(6,1);
    clinv_thick = zeros(6,1);
    clfs_thick  = zeros(6,1);
   
    %interpolate the values to the different thicknesses
    for k=1:6 % k indicate the airfoil
        clthick(k)=interp1(aoa(:,k),cl(:,k),angle_of_attack);
        cdthick(k)=interp1(aoa(:,k),cd(:,k),angle_of_attack);
        cmthick(k)=interp1(aoa(:,k),cm(:,k),angle_of_attack);
        fstat_thick(k) = interp1(aoa(:,k),f_stat(:,k),angle_of_attack);
        clinv_thick(k) = interp1(aoa(:,k),cl_inv(:,k),angle_of_attack);
        clfs_thick(k) = interp1(aoa(:,k),cl_fs(:,k),angle_of_attack);
    end

    % then interpolate to the actual thickness
    %i indicates the element nr.
    clift=interp1(thick_prof,clthick(:),thickness);
    cdrag=interp1(thick_prof,cdthick(:),thickness);
    cmom =interp1(thick_prof,cmthick(:),thickness);
    fstat =interp1(thick_prof,fstat_thick(:),thickness);
    clinv =interp1(thick_prof,clinv_thick(:),thickness);
    clfs  =interp1(thick_prof,clfs_thick(:),thickness);
end