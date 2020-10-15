function [varargout]=MyUniqueRows(varargin)

B=cell(1,nargin);
for i = 1:nargin
     B{i} = varargin{i};
end

B{1}=round(B{1},3);

[n,~]=size(B{1}); %#ok<*ASGLU>
i=1;
while i<=n-1
    j=i+1;
    while j<=n
        if B{1}(i,:)==B{1}(j,:)
            for k = 1:nargin
                 B{k}(j,:) = [];
            end
            [n,~]=size(B{1});
        else
            j=j+1;
        end
    end
    i=i+1;
end
varargout=cell(1,nargin);
for i = 1:nargin
     varargout{i}=B{i};
end
end