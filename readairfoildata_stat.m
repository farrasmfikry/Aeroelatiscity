function [aoa,cl,cd,cm,f_stat,cl_inv,cl_fs]=...
    readairfoildata_stat(file_path, ...
    filename1,filename2,...
    filename3,filename4, ...
    filename5,filename6)
strucdata=importdata([file_path filename1]);
aoa(:,1)=strucdata(:,1);
cl(:,1)=strucdata(:,2);
cd(:,1)=strucdata(:,3);
cm(:,1)=strucdata(:,4);
f_stat(:,1) =strucdata(:,5);
cl_inv(:,1) =strucdata(:,6);
cl_fs(:,1) =strucdata(:,7);

strucdata=importdata([file_path filename2]);
aoa(:,2)=strucdata(:,1);
cl(:,2)=strucdata(:,2);
cd(:,2)=strucdata(:,3);
cm(:,2)=strucdata(:,4);
f_stat(:,2) =strucdata(:,5);
cl_inv(:,2) =strucdata(:,6);
cl_fs(:,2) =strucdata(:,7);

strucdata=importdata([file_path filename3]);
aoa(:,3)=strucdata(:,1);
cl(:,3)=strucdata(:,2);
cd(:,3)=strucdata(:,3);
cm(:,3)=strucdata(:,4);
f_stat(:,3) =strucdata(:,5);
cl_inv(:,3) =strucdata(:,6);
cl_fs(:,3) =strucdata(:,7);

strucdata=importdata([file_path filename4]);
aoa(:,4)=strucdata(:,1);
cl(:,4)=strucdata(:,2);
cd(:,4)=strucdata(:,3);
cm(:,4)=strucdata(:,4);
f_stat(:,4) =strucdata(:,5);
cl_inv(:,4) =strucdata(:,6);
cl_fs(:,4) =strucdata(:,7);

strucdata=importdata([file_path filename5]);
aoa(:,5)=strucdata(:,1);
cl(:,5)=strucdata(:,2);
cd(:,5)=strucdata(:,3);
cm(:,5)=strucdata(:,4);
f_stat(:,5) =strucdata(:,5);
cl_inv(:,5) =strucdata(:,6);
cl_fs(:,5) =strucdata(:,7);

strucdata=importdata([file_path filename6]);
aoa(:,6)=strucdata(:,1);
cl(:,6)=strucdata(:,2);
cd(:,6)=strucdata(:,3);
cm(:,6)=strucdata(:,4);
f_stat(:,6) =strucdata(:,5);
cl_inv(:,6) =strucdata(:,6);
cl_fs(:,6) =strucdata(:,7);
     end