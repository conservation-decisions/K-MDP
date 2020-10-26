% Author: iadine.Chades@csiro.au
% input: stateid number
% output: site variable corresponding to the configuration of the J sites
% note: we use a ternary base.

function site=getSite(stateid,J)   
    baseTern=(3*ones(1,J)).^[0:J-1];
    site=zeros(1,J);
    for i=J:-1:1
        if stateid-2*baseTern(i)>=0
            site(i)=2;
            stateid=stateid-2*baseTern(i);
        elseif stateid-baseTern(i)>=0
            site(i)=1;
            stateid=stateid-baseTern(i);
        end
    end
end