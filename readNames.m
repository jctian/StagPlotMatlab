% read output fields and total step

function [fields, min_frame, max_frame] = readNames(directory)

cd(directory)
% in a +op folder now

clear 
listing = dir();

leng = size(listing,1);
filenames = {listing.name};
max_frame = 0;
min_frame = 10000000;
fields = [" "];

for i = 1:leng
   name = char(filenames(i));
   if contains(name, '.')
       continue
   else
       newname = split(name, '_');
       newname = string(newname(2));
       field = split(newname, {'1','2','3','4','5','6','7','8','9','0'}); % remove digites, get field
       field = string(field(1));
       if contains(fields, field) == false
           fields = fields + field + " ";
       end
       
       
       frame = split(newname, field); % {0*0cell, frame}, so second element is number of frame
       frame = str2double(frame(2));
       if(frame>max_frame)
           max_frame = frame;
       end
       if(frame<min_frame)
           min_frame = frame;
       end
   end
       
    
end




