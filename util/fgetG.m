function L=fgetG(F,n)
% Runs Fgets but rips off the carriage return and adds a space.
% Returns -1 if at the end of the file
% 
% If given an argument n, fgetg will pad L with spaces so length(L)>=n.


L=fgets(F);  % if L~=-1,     L=[L(1:end-1),' ']; end;

  if L~=-1,
      Ldec=unicode2native(L);
      index=find(Ldec==13);
      if index
          L=[L(1:index-1),' '];
      else
          L=[L(1:end-1),' '];
      end

	% Make sure that the minimum length of L is n.
	if nargin==2
   	   if length(L)< n,
           L=[L,blanks(n-length(L))];
       end
    end

  end



%alternatively, but not better

   %if nargin~=2, n=1, 
   %else          n=max(0,n-length(L));
   %end
   %L=[L,blanks(n-length(L))];  