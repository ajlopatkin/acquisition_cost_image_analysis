function maxGR = growth_manual(time, timecourse, xminus, xplus)
 
midTC = max(timecourse)./2;
ind = find(timecourse >=midTC ,1 ,'first');

if ind <= xminus
    ind = xminus+1;
end
maxGR = (timecourse(ind+xplus) - timecourse(ind-xminus)) ./ (time(ind+xplus) - time(ind-xminus));
% plot(linspace(time(ind-xminus),time(ind+xplus),20),linspace(timecourse(ind-xminus),timecourse(ind+xplus),20),'k','linewidth',4.0)

end

