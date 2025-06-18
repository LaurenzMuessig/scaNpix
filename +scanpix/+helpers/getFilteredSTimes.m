function st_filt = getFilteredSTimes(st,validPos,Fs)

%%
arguments
    st {mustBeNumeric} 
    validPos {mustBeNumericOrLogical} 
    Fs (1,1) {mustBeNumeric} = 50 %
end

%%

st_binned = ceil(st .* Fs);
ind       = validPos(st_binned);
st_filt   = st(ind);


end