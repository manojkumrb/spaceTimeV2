function observed_factor = observed_factor_setup(name_str,datain)

% Extract data for use as observed factors (Note, these cannot include NaN during estimation period
for i = 1:size(name_str);
    ustr = char(name_str(i));
    j = colnumber(ustr,datain.bpnamevec);
    if i == 1;
        w = datain.bpdata(:,j);
    else;
        w = [w datain.bpdata(:,j)];
    end;
end;
observed_factor.fac = w;
observed_factor.names = name_str;

end