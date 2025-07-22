function W = dcf_pipe(F,niter)
% returns density compensation weighting Fatrix for non-uniform nufft F
% using Jim Pipe's method
% reference:
% Pipe JG, Menon P. Sampling density compensation in MRI: rationale and an
% iterative numerical solution. Magn Reson Med. 1999 Jan;41(1):179-86.
% doi: 10.1002/(sici)1522-2594(199901)41:1<179::aid-mrm25>3.0.co;2-v.
% by David Frey (djfrey@umich.edu)
%
% inputs:
% F - NUFFT operator
% niter - number of iterations
%
% outputs:
% W - density compensation weighting operator (diagonal fatrix)
%

    if nargin < 2 || isempty(niter)
        niter = 3;
    end
    wi =  ones(size(F,1),1);
    for itr = 1:niter
        wd = real( F.arg.st.interp_table(F.arg.st, ...
            F.arg.st.interp_table_adj(F.arg.st, wi) ) );
        wi = wi ./ wd;
    end
    W = Gdiag(wi / sum(abs(wi)));
    
end