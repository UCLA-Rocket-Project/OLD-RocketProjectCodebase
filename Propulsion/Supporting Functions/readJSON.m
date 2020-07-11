fid = fopen("./Configurations/config1",'r');
config = jsondecode(fscanf(fid,'%c'));
fclose(fid);