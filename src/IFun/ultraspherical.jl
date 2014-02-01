export ultraconversion!,ultraint!

## Start of support for UFun

# diff from T -> U
ultradiff(v::Vector)=[1:length(v)-1].*v[2:end]

#int from U ->T
ultraint(v::Vector)=[0,v./[1:length(v)]]

#TODO: what about missing truncation?
function ultraint!(v::Array{Float64,2})
    for k=size(v,1):-1:2
        for j=1:size(v,2)
            @inbounds v[k,j] = v[k-1,j]/(k-1)
        end
    end
    
    for j=1:size(v)[2]
        @inbounds v[1,j] = 0.
    end
    
    v
end

function ultraint!(v::Vector{Float64})
    for k=length(v):-1:2
        @inbounds v[k] = v[k-1]/(k-1)
    end
    
    @inbounds v[1] = 0.
    
    v
end


# Convert from U -> T
function ultraiconversion(v::Vector{Float64})
    n = length(v)
    w = Array(Float64,n)
        
    if n == 1
        w[1] = v[1]
    elseif n == 2
        w[1] = v[1]
        w[2] = 2v[2]
    else
        w[end] = 2v[end];
        w[end-1] = 2v[end-1];
        
        for k = n-2:-1:2
            w[k] = 2*(v[k] + .5w[k+2]);
        end
        
        w[1] = v[1] + .5w[3];
    end
    
    w
end


# Convert T -> U
function ultraconversion(v::Vector{Float64})
    n = length(v);
    w = Array(Float64,n);
    
    if n == 1
        w[1] = v[1];
    elseif n == 2
        w[1] = v[1];
        w[2] = .5v[2];
    else
        w[1] = v[1] - .5v[3];        
    
        w[2:n-2] = .5*(v[2:n-2] - v[4:n]);
    
        w[n-1] = .5v[n-1];
        w[n] = .5v[n];        
    end
    
    w
end

function ultraconversion!(v::Vector{Float64})
    n = length(v) #number of coefficients

    if n == 1
        #do nothing
    elseif n == 2
        @inbounds v[2] *= .5
    else
        @inbounds v[1] -= .5v[3];        
    
        for j=2:n-2
            @inbounds v[j] = .5*(v[j] - v[j+2]);
        end
        @inbounds v[n-1] *= .5;
        @inbounds v[n] *= .5;                
    end
    
    return v
end

function ultraconversion!(v::Array{Float64,2})
    n = size(v)[1] #number of coefficients
    m = size(v)[2] #number of funs


    if n == 1
        #do nothing
    elseif n == 2
        for k=1:m
            @inbounds v[2,k] *= .5
        end
    else
        for k=1:m
            @inbounds v[1,k] -= .5v[3,k];        
        
            for j=2:n-2
                @inbounds v[j,k] = .5*(v[j,k] - v[j+2,k]);
            end
            @inbounds v[n-1,k] *= .5;
            v[n,k] *= .5;                
        end
    end
    
    return v
end

ultraiconversion(v::Vector{Complex{Float64}})=ultraiconversion(real(v)) + ultraiconversion(imag(v))*1.im
ultraconversion(v::Vector{Complex{Float64}})=ultraconversion(real(v)) + ultraconversion(imag(v))*1.im

