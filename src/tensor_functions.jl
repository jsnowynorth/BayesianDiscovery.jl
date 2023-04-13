


"""
tensor_mult(G, A, B, C)

Used to multiply a tensor G ∈ R(I1, I2, I3) with matrices A ∈ R(P, I1), B ∈ R(Q, I2), C ∈ R(R, I3)

"""
function tensor_mult(G, A, B, C)
  # P = size(A, 2)
  # Q = size(B, 2)
  R = size(C, 2)

  tmp = cat(dims = [3], [(A * G[:,:,r]) for r in 1:R]...)
  tmp = cat(dims = [3], [tmp[:,:,r] * B' for r in 1:R]...)

  if size(A,1) < size(B,1)
    tmp = permutedims(cat(dims = [3], [tmp[i,:,:] * C' for i in 1:size(A,1)]...), [3,1,2])
  else
    tmp = permutedims(cat(dims = [3], [tmp[:,j,:] * C' for j in 1:size(B,1)]...), [1,3,2])
  end

  return tmp
  
end


"""
unfold3(Y)

Used to get the mode-3 matrix of the tensor Y.

"""
function unfold3(Y::Array{Float64,3})

  I = size(Y, 1)
  J = size(Y, 2)
  K = size(Y, 3)

  outprod = Array{Float64}(undef, K, I * J)
  for i in 1:K
    outprod[i, :] = copy(reshape(Y[:, :, i], :)')
  end

  return outprod

end

function unfold3(Y::Array{Union{Missing,Float64},3})

  I = size(Y, 1)
  J = size(Y, 2)
  K = size(Y, 3)

  outprod = Array{Union{Missing,Float64}}(undef, K, I * J)
  for i in 1:K
    outprod[i, :] = copy(reshape(Y[:, :, i], :)')
  end

  return outprod

end


"""
fold3(Y)

Used to construct a tensor from the mode-3 matrix Y.

## Examples

Y3 = unfold3(Y)
I = size(Y, 1)
J = size(Y, 2)
K = size(Y, 3)

Y = fold3(Y3, I, J, K)

"""
function fold3(Y, I, J, K)

    outprod = Array{Float64}(undef, I, J, K)
    for i in 1:K
        outprod[:,:,i] = reshape(Y[i,:], I, J)
    end

    return outprod
    
end