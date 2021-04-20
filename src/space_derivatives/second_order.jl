###############################
#
#        SECOND ORDER
#
###############################

# Thanks to https://web.media.mit.edu/~crtaylor/calculator.html
# For different stencil sizes:
# c1 = [1, -2, 1]
# c2 = [-1, 16, -30, 16, -1] / 12
# c3 = [2, -27, 270, -490, 270, -27, 2] / 180
# c4 = [-9, 128, -1008, 8064, -14350, 8064, …] / 5040
# etc

function fd_diff2_ord2!(dx2_U, U, dx⁻¹)
    local c1 = SA[1.0, -2.0, 1.0]
    for i = 2 : length(U)-1
        dx2_U[i] = dot(c1, U[i-1 : i+1])
    end
    dx2_U[1] = dx2_U[2]
    dx2_U[end] = dx2_U[end-1]
    
    @. dx2_U *= dx⁻¹ ^ 2

    nothing
end

function fd_diff2_ord4!(dx2_U, U, dx⁻¹)
    local c1 = SA[1.0, -2.0, 1.0]
    local c2 = SA[-1.0, 16.0, -30.0, 16.0, -1.0] * 0.08333333333333333333
    for i in 3 : length(U)-2
        dx2_U[i] = dot(c2, U[i-2 : i+2])
    end
    
    dx2_U[2] = dot(c1, U[1:3])
    dx2_U[end-1] = dot(c1, U[end-2 : end])
    
    dx2_U[1] = dx2_U[2]
    dx2_U[end] = dx2_U[end-1]

    @. dx2_U *= dx⁻¹ ^ 2

    nothing
end

function fd_diff2_ord6!(dx2_U, U, dx⁻¹)
    local c1 = SA[1.0, -2.0, 1.0]
    local c2 = SA[-1.0, 16.0, -30.0, 16.0, -1.0] * 0.08333333333333333333
    local c3 = SA[2.0, -27.0, 270.0, -490.0, 270.0, -27.0, 2.0]
    c3 *= 0.00555555555555555555
    for i in 4 : length(U)-3
        dx2_U[i] = dot(c3, U[i-3 : i+3])
    end

    dx2_U[3] = dot(c2, U[1:5])
    dx2_U[end-2] = dot(c2, U[end-4 : end])
    
    dx2_U[2] = dot(c1, U[1:3])
    dx2_U[end-1] = dot(c1, U[end-2 : end])
    
    dx2_U[1] = 2.0*dx2_U[2] - dx2_U[3]
    dx2_U[end] = 2.0*dx2_U[end-1] - dx2_U[end-2]

    @. dx2_U *= dx⁻¹ ^ 2

    nothing
end
