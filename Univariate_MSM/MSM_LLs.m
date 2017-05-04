
function LLs= MSM_LLs(input,kbar,n,data,A_template)


if input(3) >=1
    input(3)=0.9999;
end

if input(2) >=2
    input(2)=1.9999;
end

[~,LLs] = MSM_likelihood(input,kbar,n, data,A_template);
  

end
