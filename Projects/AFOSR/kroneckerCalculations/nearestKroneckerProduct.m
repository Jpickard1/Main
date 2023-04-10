function [ B, C, D ] = nearestKroneckerProduct( A, SizeB, SizeC, Hermitian ) %#codegen
% REFERENCE: https://doi.org/10.1007/978-94-015-8196-7_17
%
% Auth: Tom Holden
%       thomas.holden@gmail.com
% Date: June 22, 2018
%
% ORIGINIAL CODE: https://gist.github.com/tholden/58dd9a8991daa750ae36a633fe7060a4
    m = size( A, 1 );
    n = size( A, 2 );
    m1 = SizeB( 1 );
    n1 = SizeB( 2 );
    m2 = SizeC( 1 );
    n2 = SizeC( 2 );
    
    if nargin < 4
        Hermitian = false;
    end

    assert( m == m1 * m2 );
    assert( n == n1 * n2 );
    
    if Hermitian
        assert( m1 == n1 );
        assert( m2 == n2 );
        A = 0.5 * ( A + A' );
    end
    
    R = reshape( permute( reshape( A, [ m2, m1, n2, n1 ] ), [ 2 4 1 3 ] ), m1 * n1, m2 * n2 );
    
    [ B, S, C ] = svds( R, 1 );
    
    SqrtS = sqrt( S );
    
    B = reshape( B * SqrtS, m1, n1 );
    C = reshape( C * SqrtS, m2, n2 );
    
    if Hermitian
        B = 0.5 * ( B + B' );
        C = 0.5 * ( C + C' );
        
        if all( diag( B ) < 0 ) && all( diag( C ) < 0 )
            B = -B;
            C = -C;
        end
    end
    
    if nargout > 2
        D = A - kron( B, C );
        if Hermitian
            D = 0.5 * ( D + D' );
        end
    end
end