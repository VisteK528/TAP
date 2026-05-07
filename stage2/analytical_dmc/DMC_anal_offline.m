function [ Ke , Ku , Kz ] = DMC_anal_offline (S , S_z , D , D_z , N , Nu , lambda , psi )
    [~ , ny , nu ] = size ( S ) ;
    M = DMC_M (S ,N , Nu ) ;
    Mp = DMC_Mp (S ,D , N ) ;
    
    if D_z > 0
        Mpz = DMC_Mpz ( S_z , D_z , N ) ;
    end
    
    Psi = kron ( eye ( N ) , diag ( psi ) ) ;
    Lambda = kron ( eye ( Nu ) , diag ( lambda ) ) ;
    K =( M'* Psi * M + Lambda ) \M'* Psi ;
    K1 = K (1: nu ,:) ;
    Ke = 0;
    
    for i = 1: N
        Ke = Ke + K1 (: ,1+( i -1) * ny : i * ny ) ;
    end
    
    Ku = K1 * Mp ;
    
    if D_z > 0
        Kz = K1 * Mpz ;
    else
        Kz = [];
    end
end