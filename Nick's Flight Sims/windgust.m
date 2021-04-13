function vgust = windgust(gust,t,tgust,dtgust,dtrise)
    if t >= tgust && t <= tgust+dtgust+2*dtrise
        if t <= tgust+dtrise
            if dtrise == 0
                vgust = gust;
            else
                vgust = gust/2*(1-cos(pi*(t-tgust)/dtrise));
            end
        elseif t <= tgust+dtrise+dtgust
            vgust = gust;
        else
            if dtrise == 0
                vgust = gust;
            else
                vgust = gust/2*(1-cos(pi*(t-tgust-dtgust)/dtrise));
            end
        end    
    else
        vgust = 0;
    end
end