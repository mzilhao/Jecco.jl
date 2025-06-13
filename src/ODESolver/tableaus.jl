struct RK2Tableau{T} <: Tableau
    a21 :: T
    b1  :: T
    b2  :: T
    c2  :: T
end

struct RK4Tableau{T} <: Tableau
    a21 :: T
    a31 :: T
    a32 :: T
    a41 :: T
    a42 :: T
    a43 :: T
    b1  :: T
    b2  :: T
    b3  :: T
    b4  :: T
    c2  :: T
    c3  :: T
    c4  :: T
end

function RK2Tableau(T::Type = Float64)
    a21 = convert(T, 1//2)
    b1  = convert(T, 0)
    b2  = convert(T, 1)
    c2  = convert(T, 1//2)
    RK2Tableau{T}(a21, b1, b2, c2)
end

function RK4Tableau(T::Type = Float64)
    a21 = convert(T, 1//2)
    a31 = convert(T, 0)
    a32 = convert(T, 1//2)
    a41 = convert(T, 0)
    a42 = convert(T, 0)
    a43 = convert(T, 1)
    b1  = convert(T, 1//6)
    b2  = convert(T, 1//3)
    b3  = convert(T, 1//3)
    b4  = convert(T, 1//6)
    c2  = convert(T, 1//2)
    c3  = convert(T, 1//2)
    c4  = convert(T, 1)
    RK4Tableau{T}(a21, a31, a32, a41, a42, a43, b1, b2, b3, b4, c2, c3, c4)
end
