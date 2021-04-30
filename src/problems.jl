using OffsetArrays


Flux = OffsetVector
zero_flux(Nx :: Int) = Flux(zeros(Nx+1), 0 : Nx)
zero_flux(U :: Vector) = Flux(zeros(length(U)+1), 0 : length(U))

struct ScalarQuantity
    data :: Vector
    posFlux :: Flux
    negFlux :: Flux

    ScalarQuantity(U :: Vector) = new(copy(U), zero_flux(U), zero_flux(U))
end


struct RelaxPb
    U :: ScalarQuantity
    V :: ScalarQuantity

    fluxP :: Flux
    fluxV :: Flux
    dxP :: Vector
    dxV :: Vector

    f :: Function
    p :: Function

    function RelaxPb(U :: Vector, V :: Vector,
                     f :: Function, p :: Function)

        Nx = length(U)
        pbU, pbV = ScalarQuantity(U), ScalarQuantity(V)
        fluxV, fluxP = zero_flux(Nx), zero_flux(Nx)
        dxV, dxP = zeros(Nx), zeros(Nx)

        new(pbU,pbV,fluxP,fluxV,dxV,dxP,f,p)
    end
end

function Base.copy(pb :: RelaxPb)
    RelaxPb(pb.U.data,pb.V.data,pb.f,pb.p)
end