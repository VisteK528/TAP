function [ dU ]= DMC_anal_dU (y , y_zad , Ke , Ku , Kz , dUp , dZ )
    dU = Ke *( y_zad - y ) ;
    nu = size ( Ke ,1) ;
    
    for i = 1: size ( Ku ,2) / nu
        dU = dU - Ku (: ,1+( i -1) * nu : i * nu ) * dUp (: , i ) ;
    end
    
    if ~ isempty ( Kz )
        for i = 1: size ( Kz ,2)
            dU = dU - Kz (: , i ) * dZ (: , i ) ;
        end
    end
end