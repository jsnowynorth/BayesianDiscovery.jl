
# U through U³
function dU(A, ϕ, ψ)
    (ϕ ⊗ ψ)'
end

function dU²(A, ϕ, ψ)
    2 .* (A * (ϕ ⊗ ψ)) .*(ϕ ⊗ ψ)'
end

function dU³(A, ϕ, ψ)
    3 .* (A * (ϕ ⊗ ψ)).^2 .* (ϕ ⊗ ψ)'
end


# U_x through U_xyy
function dU_x(A, ϕ, ψ_x)
    (ϕ ⊗ ψ_x)'
end

function dU_y(A, ϕ, ψ_y)
    (ϕ ⊗ ψ_y)'
end

function dU_xx(A, ϕ, ψ_xx)
    (ϕ ⊗ ψ_xx)'
end

function dU_yy(A, ϕ, ψ_yy)
    (ϕ ⊗ ψ_yy)'
end

function dU_xy(A, ϕ, ψ_xy)
    (ϕ ⊗ ψ_xy)'
end

function dU_xxx(A, ϕ, ψ_xxx)
    (ϕ ⊗ ψ_xxx)'
end

function dU_yyy(A, ϕ, ψ_yyy)
    (ϕ ⊗ ψ_yyy)'
end

function dU_xxy(A, ϕ, ψ_xxy)
    (ϕ ⊗ ψ_xxy)'
end

function dU_xyy(A, ϕ, ψ_xyy)
    (ϕ ⊗ ψ_xyy)'
end


# UU_x through UU_xyy
function dUU_x(A, ϕ, ψ, ψ_x)
    (((ϕ ⊗ ψ)') .* (A * (ϕ ⊗ ψ_x))) + ((A * (ϕ ⊗ ψ)) .* ((ϕ ⊗ ψ_x)'))
end

function dUU_y(A, ϕ, ψ, ψ_y)
    (((ϕ ⊗ ψ)') .* (A * (ϕ ⊗ ψ_y))) + ((A * (ϕ ⊗ ψ)) .* ((ϕ ⊗ ψ_y)'))
end

function dUU_xx(A, ϕ, ψ, ψ_xx)
    (((ϕ ⊗ ψ)') .* (A * (ϕ ⊗ ψ_xx))) + ((A * (ϕ ⊗ ψ)) .* ((ϕ ⊗ ψ_xx)'))
end

function dUU_yy(A, ϕ, ψ, ψ_yy)
    (((ϕ ⊗ ψ)') .* (A * (ϕ ⊗ ψ_yy))) + ((A * (ϕ ⊗ ψ)) .* ((ϕ ⊗ ψ_yy)'))
end

function dUU_xy(A, ϕ, ψ, ψ_xy)
    (((ϕ ⊗ ψ)') .* (A * (ϕ ⊗ ψ_xy))) + ((A * (ϕ ⊗ ψ)) .* ((ϕ ⊗ ψ_xy)'))
end

function dUU_xxx(A, ϕ, ψ, ψ_xxx)
    (((ϕ ⊗ ψ)') .* (A * (ϕ ⊗ ψ_xxx))) + ((A * (ϕ ⊗ ψ)) .* ((ϕ ⊗ ψ_xxx)'))
end

function dUU_yyy(A, ϕ, ψ, ψ_yyy)
    (((ϕ ⊗ ψ)') .* (A * (ϕ ⊗ ψ_yyy))) + ((A * (ϕ ⊗ ψ)) .* ((ϕ ⊗ ψ_yyy)'))
end

function dUU_xxy(A, ϕ, ψ, ψ_xxy)
    (((ϕ ⊗ ψ)') .* (A * (ϕ ⊗ ψ_xxy))) + ((A * (ϕ ⊗ ψ)) .* ((ϕ ⊗ ψ_xxy)'))
end

function dUU_xyy(A, ϕ, ψ, ψ_xyy)
    (((ϕ ⊗ ψ)') .* (A * (ϕ ⊗ ψ_xyy))) + ((A * (ϕ ⊗ ψ)) .* ((ϕ ⊗ ψ_xyy)'))
end


# U²U_x through U²U_xyy
function dU²U_x(A, ϕ, ψ, ψ_x)
    ((2 .* A * (ϕ ⊗ ψ) * (ϕ ⊗ ψ)') .* (A * (ϕ ⊗ ψ_x))) + (((A * (ϕ ⊗ ψ)) .^2) .* (ϕ ⊗ ψ_x)')
end

function dU²U_y(A, ϕ, ψ, ψ_y)
    ((2 .* A * (ϕ ⊗ ψ) *(ϕ ⊗ ψ)') .* (A * (ϕ ⊗ ψ_y))) + (((A * (ϕ ⊗ ψ)) .^2) .* (ϕ ⊗ ψ_y)')
end

function dU²U_xx(A, ϕ, ψ, ψ_xx)
    ((2 .* A * (ϕ ⊗ ψ) * (ϕ ⊗ ψ)') .* (A * (ϕ ⊗ ψ_xx))) + (((A * (ϕ ⊗ ψ)) .^2) .* (ϕ ⊗ ψ_xx)')
end

function dU²U_yy(A, ϕ, ψ, ψ_yy)
    ((2 .* A * (ϕ ⊗ ψ) * (ϕ ⊗ ψ)') .* (A * (ϕ ⊗ ψ_yy))) + (((A * (ϕ ⊗ ψ)) .^2) .* (ϕ ⊗ ψ_yy)')
end

function dU²U_xy(A, ϕ, ψ, ψ_xy)
    ((2 .* A * (ϕ ⊗ ψ) * (ϕ ⊗ ψ)') .* (A * (ϕ ⊗ ψ_xy))) + (((A * (ϕ ⊗ ψ)) .^2) .* (ϕ ⊗ ψ_xy)')
end

function dU²U_xxx(A, ϕ, ψ, ψ_xxx)
    ((2 .* A * (ϕ ⊗ ψ) * (ϕ ⊗ ψ)') .* (A * (ϕ ⊗ ψ_xxx))) + (((A * (ϕ ⊗ ψ)) .^2) .* (ϕ ⊗ ψ_xxx)')
end

function dU²U_yyy(A, ϕ, ψ, ψ_yyy)
    ((2 .* A * (ϕ ⊗ ψ) * (ϕ ⊗ ψ)') .* (A * (ϕ ⊗ ψ_yyy))) + (((A * (ϕ ⊗ ψ)) .^2) .* (ϕ ⊗ ψ_yyy)')
end

function dU²U_xxy(A, ϕ, ψ, ψ_xxy)
    ((2 .* A * (ϕ ⊗ ψ) * (ϕ ⊗ ψ)') .* (A * (ϕ ⊗ ψ_xxy))) + (((A * (ϕ ⊗ ψ)) .^2) .* (ϕ ⊗ ψ_xxy)')
end

function dU²U_xyy(A, ϕ, ψ, ψ_xyy)
    ((2 .* A * (ϕ ⊗ ψ) * (ϕ ⊗ ψ)') .* (A * (ϕ ⊗ ψ_xyy))) + (((A * (ϕ ⊗ ψ)) .^2) .* (ϕ ⊗ ψ_xyy)')
end


# U³U_x through U³U_xyy
function dU³U_x(A, ϕ, ψ, ψ_x)
    (3 .*(A * (ϕ ⊗ ψ)) .^2 * (ϕ ⊗ ψ)') .* (A * (ϕ ⊗ ψ_x)) + (((A * (ϕ ⊗ ψ)) .^3) .* (ϕ ⊗ ψ_x)')
end

function dU³U_y(A, ϕ, ψ, ψ_y)
    (3 .*(A * (ϕ ⊗ ψ)) .^2 * (ϕ ⊗ ψ)') .* (A * (ϕ ⊗ ψ_y)) + (((A * (ϕ ⊗ ψ)) .^3) .*(ϕ ⊗ ψ_y)')
end

function dU³U_xx(A, ϕ, ψ, ψ_xx)
    (3 .*(A * (ϕ ⊗ ψ)) .^2 * (ϕ ⊗ ψ)') .* (A * (ϕ ⊗ ψ_xx)) + (((A * (ϕ ⊗ ψ)) .^3) .* (ϕ ⊗ ψ_xx)')
end

function dU³U_yy(A, ϕ, ψ, ψ_yy)
    (3 .*(A * (ϕ ⊗ ψ)) .^2 * (ϕ ⊗ ψ)') .* (A * (ϕ ⊗ ψ_yy)) + (((A * (ϕ ⊗ ψ)) .^3) .* (ϕ ⊗ ψ_yy)')
end

function dU³U_xy(A, ϕ, ψ, ψ_xy)
    (3 .*(A * (ϕ ⊗ ψ)) .^2 * (ϕ ⊗ ψ)') .* (A * (ϕ ⊗ ψ_xy)) + (((A * (ϕ ⊗ ψ)) .^3) .*(ϕ ⊗ ψ_xy)')
end

function dU³U_xxx(A, ϕ, ψ, ψ_xxx)
    (3 .*(A * (ϕ ⊗ ψ)) .^2 * (ϕ ⊗ ψ)') .* (A * (ϕ ⊗ ψ_xxx)) + (((A * (ϕ ⊗ ψ)) .^3) .* (ϕ ⊗ ψ_xxx)')
end

function dU³U_yyy(A, ϕ, ψ, ψ_yyy)
    (3 .*(A * (ϕ ⊗ ψ)) .^2 *(ϕ ⊗ ψ)') .* (A * (ϕ ⊗ ψ_yyy)) + (((A * (ϕ ⊗ ψ)) .^3) .* (ϕ ⊗ ψ_yyy)')
end

function dU³U_xxy(A, ϕ, ψ, ψ_xxy)
    (3 .*(A * (ϕ ⊗ ψ)) .^2 * (ϕ ⊗ ψ)') .* (A * (ϕ ⊗ ψ_xxy)) + (((A * (ϕ ⊗ ψ)) .^3) .* (ϕ ⊗ ψ_xxy)')
end

function dU³U_xyy(A, ϕ, ψ, ψ_xyy)
    (3 .*(A * (ϕ ⊗ ψ)) .^2 * (ϕ ⊗ ψ)') .* (A * (ϕ ⊗ ψ_xyy)) + (((A * (ϕ ⊗ ψ)) .^3) .* (ϕ ⊗ ψ_xyy)')
end

# UxU_xx through UxxU_xxx
function dUxU_xx(A, ϕ, ψ_x, ψ_xx)
    (ϕ ⊗ ψ_x)' .* (A * (ϕ ⊗ ψ_xx)) + (A * (ϕ ⊗ ψ_x)) .* (ϕ ⊗ ψ_xx)'
end

function dUxU_xxx(A, ϕ, ψ_x, ψ_xxx)
    (ϕ ⊗ ψ_x)' .* (A * (ϕ ⊗ ψ_xxx)) + (A * (ϕ ⊗ ψ_x)) .* (ϕ ⊗ ψ_xxx)'
end

function dUxxU_xxx(A, ϕ, ψ_xx, ψ_xxx)
    (ϕ ⊗ ψ_xx)' .* (A * (ϕ ⊗ ψ_xxx)) + (A * (ϕ ⊗ ψ_xx)) .* (ϕ ⊗ ψ_xxx)'
end

# U * X through U_xxx * X
function dUX(A, ϕ, ψ, X)
    (ϕ ⊗ ψ)' .* X
end

function dU_xX(A, ϕ, ψ_x, X)
    (ϕ ⊗ ψ_x)' .* X
end

function dU_xxX(A, ϕ, ψ_xx, X)
    (ϕ ⊗ ψ_xx)' .* X
end

function dU_xxxX(A, ϕ, ψ_xxx, X)
    (ϕ ⊗ ψ_xxx)' .* X
end