function [dates,slcnames,dn,nd]=list_slcs(params)

slcdir=params.slcdir;

files=dir([slcdir '2*']);
dates=[];
for i=1:length(files)
    dates(i).name=files(i).name(1:8);
    dates(i).dn=datenum(dates(i).name,'yyyymmdd');
    slcnames{i}=[slcdir dates(i).name '/' dates(i).name '.slc.full'];
 end

