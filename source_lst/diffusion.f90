MODULE properties 

!****** Description: 
!****** This model is using to 
!****** calculate the diffusivity properties(mi_k_D)  of the 
!****** different species in function of  
!****** the temperature. 
!****** Author: Jhonatan 
!****** Date: 16/11/2017   
!paper:Transport Coefficients for the NASA Çewis Chemical Equilibrium Program
!author: Roger Svehla 1995

use global

contains

!subroutine molecular_properties()
!implicit none
!integer :: i,j
!real(kind=ip)::	a_mi,b_mi,c_mi,d_mi
!real(kind=ip)::	a_k,b_k,c_k,d_k
!real(kind=ip)::	a_cp,b_cp,c_cp,d_cp,e_cp
!real(kind=ip),dimension(i_s-g_p:d_s+g_p,s_s-g_p:e_s+g_p)::T_ins
!
!
!call  visc_N(a_mi,b_mi,c_mi,d_mi)
!call  cond_N(a_k ,b_k ,c_k ,d_k )
!call  heatcapacity_N(a_cp ,b_cp ,c_cp ,d_cp,e_cp,Mw1)
!
!do i=i_s,d_s
!    do j=s_s,e_s
!
!	!Temperatura instatantanea dimensional
!
!	T_ins(i,j) =	Tempba(i,j)	+	Temp(i,j)
!	   
!      mi(i,j)   = 1.d0/mi_0*(dexp(A_mi*dlog((T_ins(i,j)*T_0)) + B_mi/(T_ins(i,j)*T_0) & 
!    & 		 + C_mi/(T_ins(i,j)*T_0)/(T_ins(i,j)*T_0) + D_mi)) 
!
!      ka(i,j)   =  1.d0/k_0*(dexp(A_k*dlog((T_ins(i,j)*T_0))  + B_k/(T_ins(i,j)*T_0)  & 
!    & + C_k/(T_ins(i,j)*T_0)/(T_ins(i,j)*T_0)  + D_k)) 
!
!	Cp(i,j)	= 1.d0/cp_0*(A_cp + B_cp*(T_ins(i,j)*T_0)  + C_cp*(T_ins(i,j)*T_0)**2.d0       & 
!	&               + D_cp*(T_ins(i,j)*T_0)**3.d0 + E_cp*(T_ins(i,j)*T_0)**4.d0) 
!
!    end do  
!end do 
!END SUBROUTINE molecular_properties 
!
!subroutine init_diffusion(mib,kb,cpb,cpb1,cpb2,D_1m,Mw,Ru,Rmix,gammab,gamma_j,P_r,X1,X2,Tb,T_r,T_0,T_ref,N)
subroutine init_diffusion(mib,kb,cpb,cpb1,cpb2,D_1m,Mw,Ru,Rmix,&
                       &  gammab,gamma_ref,gamma_j,gamma_ratio,P_r,a_ref,&  
                       &  a_s,a_1,a_2,X1,X2,Tb,T_r,T_0,T_ref,N)

!IN 
!...T_0    = Characteristic temperature[k]
!...T1     = Temperature specie 1  dimensionaless
!...T2     = Temperature specie 2  dimensionaless 
!...T0     = Temperature ambient   dimensionaless 
!...Tb[j]  = Temperature base of the mixing, global.f90
!...k1     = Conductivity specie 1  dimensionaless
!...k2     = Conductivity specie 2  dimensionaless 

!OUT
!...mi_0   = Characteristic viscosity  specie 1 [gr/cm*sx10-6]
!...k_0    = Characteristic conductivity specie 1 [x10-6]
!...mib[j] = viscosity  base of the mixing, global.f90  

!Internals
!Coeficcient fitting polynome 
!paper:Transport Coefficients for the NASA Lewis Chemical Equilibrium Program
!author: Roger Svehla 1995
!A_MI
!B_MI
!C_MI
!D_MI

implicit none
integer :: i,j
real(kind=ip)::T_i,T_o,T_r,T_0,T_ref
!f2py intent(in)::T_i,T_o,T_r,T_0
integer      ::N
!f2py intent(in)::N
real(kind=ip),dimension(0:N)::Tb,X1,X2,X3,mib,kb,cpb
!f2py intent(in)::Tb,X1,X2,X3  
!f2py intent(out)::mib,kb,cpb,cpb1,cpb2 
real(kind=ip),dimension(0:N)::mib1,kb1,cpb1,mib2,kb2,cpb2,mib3,kb3,cpb3
real(kind=ip)::mi_1,k_1,cp_1,mi_2,k_2,cp_2,mi_3,k_3,cp_3
!f2py intent(out)::gamma_j,gamma_ref 
real(kind=ip)::gamma_1,gamma_2,gamma_j
real(kind=ip),dimension(0:N)::phi12_mi,phi21_mi,phi31_mi,phi13_mi,phi23_mi,phi32_mi
real(kind=ip),dimension(0:N)::phi12_k ,phi21_k ,phi31_k ,phi13_k ,phi23_k ,phi32_k
real(kind=ip),dimension(0:N)::phi12_cp,phi21_cp ,phi31_cp ,phi13_cp ,phi23_cp ,phi32_cp
real(kind=ip)::phi11,phi22,phi33,phi_1
real(kind=ip)::Mw1,Mw2,Mw3,Ru,gamma_ratio,a_1,a_2,a_ref
real(kind=ip)::mi_ref,k_ref,cp_ref,gamma_ref,Mw_ref,Rmix_ref,Rmix_1,Rmix_2 
!f2py intent(out)::Mw,Ru,Rmix,gammab,gamma_ratio,P_r,a_s,a_ref,a_1,a_2
!f2py intent(out)::Mw,Ru,Rmix,gammab,P_r
real(kind=ip),dimension(0:N)::Mw,R1,Rmix,gammab,P_r,a_s
real(kind=ip),dimension(0:N)::mib_1,mib_2,mib_3
real(kind=ip),dimension(0:N)::kb_1,kb_2,kb_3
real(kind=ip),dimension(0:N)::cpb_1,cpb_2,cpb_3
real(kind=ip),dimension(0:N)::D_1m,D_2m,D_3m
!f2py intent(out)::D_1m,D_2m,D_3m

!Molar gas constant
Ru      = 8.314510d0![J/g*K]

call proper_2(mi_1, mib1, k_1, kb1, cp_1, cpb1, Tb, T_r, T_ref, Mw1, Ru, N)
call proper_2(mi_2, mib2, k_2, kb2, cp_2, cpb2, Tb, T_0, T_ref, Mw2, Ru, N)


!write(*,*)gamma_0


phi_1   = 1.d0 
phi11   = phi_1 
phi22   = phi_1 
phi33   = phi_1 

call mix_proper(mib,X1,X2,X3,mib1,mib2,mib3,Mw1,Mw2,Mw3,N)

call mix_proper(kb ,X1,X2,X3,kb1 ,kb2 ,kb3 ,Mw1,Mw2,Mw3,N)
!


do j=s_s,e_s

        Mw(j)        =    1.d0/(X1(j)/Mw1+X2(j)/Mw2)
        Rmix(j)      =    Ru/Mw(j) 
        !Anderson: Frozem cp
        cpb (j)      =    X1(j)*cpb1(j)+X2(j)*cpb2(j)
        gammab(j)    =    1.d0/(1.d0-Rmix(j)/cpb(j))     
        P_r(j)       =    1.d0/gammab(j)
        !*1000gr
        a_s(j)       =    sqrt(gammab(j)*Rmix(j)*1000*(Tb(j)*T_ref))

end do 



call mix_proper_diffusion(D_1m,X1,X2,X3,Mw1,Mw2,Mw3,P_r,Tb,T_r,T_0,T_ref,N)

!Reference specified heat ratio 

gamma_1         =       1.d0/(1.d0-(Ru/Mw1)/cp_1)     
gamma_2         =       1.d0/(1.d0-(Ru/Mw2)/cp_2)     
Rmix_1          =        Ru/Mw1
Rmix_2          =        Ru/Mw2
a_1             =        sqrt(gamma_1*Rmix_1*1000*(T_r*T_ref))
a_2             =        sqrt(gamma_2*Rmix_2*1000*(T_0*T_ref))

mi_ref          =          mi_1
k_ref           =           k_1
cp_ref          =          cp_1
Mw_ref          =           Mw1 
a_ref           =           a_1
gamma_ref       =       gamma_1 
Rmix_ref        =     Ru/Mw_ref

!write(*,*)'a_2',a_1,'k_2',k_1,'cp_2',cp_1

gamma_ratio     = gamma_2/gamma_1 

gamma_j         = gamma_ref/(gamma_ref-1.d0)

!dimensionalized the mixture viscosity by the mi_r:

mib     =       mib/mi_ref


!dimensionalized the mixture conductivity by the k_r:
kb      =       kb/k_ref

!dimensionalized the mixture heat capacity by the cp_r:
cpb     =       cpb/cp_ref

gammab  =       gammab/gamma_ref

!dimensionalized the mixture constant:

Rmix    =       Rmix/Rmix_ref 

!non dimensional speed of sound of  the mixture 
a_s     =       a_s/a_ref 


end subroutine init_Diffusion

subroutine mix_proper_diffusion(D12,X1,X2,X3,Mw1,Mw2,Mw3,P_r,Tb,T_r,T_0,T_ref,N)

implicit none                  
integer :: i,j,N
real(kind=ip)::Mw1,Mw2,Mw3
real(kind=ip)::Mw11_df,Mw22_df,Mw33_df
real(kind=ip)::Mw12_df,Mw13_df,Mw23_df
real(kind=ip)::T_r,T_ref,T_i,T_o,T_0
real(kind=ip)::epsilonk1_df ,epsilonk2_df ,epsilonk3_df
real(kind=ip)::epsilonk12_df,epsilonk13_df,epsilonk23_df
real(kind=ip)::delta1_df,delta2_df,delta3_df
real(kind=ip)::delta11_df,delta22_df,delta33_df
real(kind=ip)::delta12_df,delta13_df,delta23_df
real(kind=ip)::A_df,B_df,C_df,D_df,E_df,F_df,G_df,H_df 
real(kind=ip),dimension(0:N)::X1,X2,X3,P_r
real(kind=ip),dimension(0:N)::omega11_df,omega22_df,omega33_df
real(kind=ip),dimension(0:N)::omega12_df,omega13_df,omega23_df
real(kind=ip),dimension(0:N)::D_1m,D_2m,D_3m
real(kind=ip),dimension(0:N)::D11,D22,D33,D12,D13,D23
real(kind=ip),dimension(0:N)::Tb,Tb_st,Tb_st1,Tb_st2,Tb_st3
real(kind=ip),dimension(0:N)::Tb_st12,Tb_st13,Tb_st23

!ALL IS MESSURE RESPECT TO T_r


call diffusion_O2(delta1_df,epsilonk1_df)
call diffusion_O2(delta2_df,epsilonk2_df)
call diffusion_O2(delta3_df,epsilonk3_df)

epsilonk12_df    =       (epsilonk1_df*epsilonk2_df)**(0.5d0)

delta12_df       =       (delta1_df+delta2_df)/(2.d0)

Mw12_df          =       2.d0*1.d0/(1.d0/Mw1+1.d0/Mw2)

A_df    =       1.06036d0
B_df    =       0.15610d0
C_df    =       0.19300d0
D_df    =       0.47635d0
E_df    =       1.03587d0
F_df    =       1.52996d0
G_df    =       1.76474d0
H_df    =       3.89411d0

Tb_st   =       Tb*T_ref 

Tb_st12 =       Tb_st/epsilonk12_df 


!Colision Integral
omega11_df      =     A_df/((Tb_st1)**(B_df)   )   &
                & +   C_df/((dexp(D_df*Tb_st1)))   &
                & +   E_df/((dexp(F_df*Tb_st1)))   &
                & +   G_df/((dexp(H_df*Tb_st1)))   

omega22_df      =     A_df/((Tb_st2)**(B_df)   )   &
                & +   C_df/((dexp(D_df*Tb_st2)))   &
                & +   E_df/((dexp(F_df*Tb_st2)))   &
                & +   G_df/((dexp(H_df*Tb_st2)))   

omega12_df      =     A_df/((Tb_st12)**(B_df)  )    &
                & +   C_df/((dexp(D_df*Tb_st12)))   &
                & +   E_df/((dexp(F_df*Tb_st12)))   &
                & +   G_df/((dexp(H_df*Tb_st12)))   

do j=s_s,e_s


        D12(j)    = 0.00266d0*Tb_st12(j)**(1.5d0)/(P_r(j)*Mw12_df**(0.5d0)*delta12_df**2.d0*omega12_df(j))


!        write(*,*)D_1m(j),X1(j),X2(j)
end do 

end subroutine mix_proper_diffusion 

subroutine mix_proper(Pmix,X1,X2,X3,Pb1,Pb2,Pb3,Mw1,Mw2,Mw3,N)

implicit none
integer :: i,j,N
real(kind=ip)::Mw1,Mw2,Mw3
real(kind=ip)::phi11,phi22,phi33
real(kind=ip),dimension(0:N)::phi12,phi13,phi21,phi23,phi31,phi32
real(kind=ip),dimension(0:N)::X1,X2,X3,Pb1,Pb2,Pb3,Pmix,Pb_1,Pb_2,Pb_3

phi11   = 1.d0 
phi22   = 1.d0 

call phi_f(phi12 ,Mw1,Mw2,Pb1,Pb2,N)
call phi_f(phi21 ,Mw2,Mw1,Pb2,Pb1,N)

do j=s_s,e_s

        Pb_1(j) = X1(j)*Pb1(j)/(X1(j)*phi11   +X2(j)*phi12(j))
        Pb_2(j) = X2(j)*Pb2(j)/(X1(j)*phi21(j)+X2(j)*phi22   )

        Pmix(j) = Pb_1(j)+Pb_2(j)
end do 


end subroutine mix_proper 

subroutine phi_f(phi_ij,Mwi,Mwj,Pbi,Pbj,N)

implicit none

integer::j,N
real(kind=ip)::Mwi,Mwj,phi_1
real(kind=ip),dimension(0:N)::phi_ij,Pbi,Pbj

phi_1   = 1.d0/dsqrt(8.d0) 

do j=s_s,e_s

        phi_ij(j)  = phi_1*(1.d0+Mwi/Mwj)**(-0.5d0)*(1.d0+(Pbi(j)/Pbj(j))**(0.5d0)*(Mwj/Mwi)**(0.25d0))**2.d0 

end do 

end subroutine phi_f 

subroutine proper_1(mi_d,mib,k_d,kb,cp_d,cpb,Tb,T_ref,T_r,Mw,Ru,N)

implicit none
integer :: i,j
real(kind=ip)::T_b,T_ref,T_r,T_d
real(kind=ip)::A_mi,B_mi,C_mi,D_mi
real(kind=ip)::A_k , B_k ,C_k ,D_k 
real(kind=ip)::A_cp ,B_cp ,C_cp ,D_cp,E_cp
integer      ::N
real(kind=ip),dimension(0:N)::Tb
real(kind=ip),dimension(0:N)::mib,kb,cpb
real(kind=ip)::mi_d,k_d,cp_d
real(kind=ip)::Mw,Ru


CALL          visc_H2(A_mi, B_mi , C_mi ,D_mi)
CALL          cond_H2(A_k , B_k  , C_k  ,D_k )
CALL  heatcapacity_H2(A_cp, B_cp , C_cp ,D_cp,E_cp,Mw)

!.....Dimensionaless equation for 
!.....temperature, using T_r T'=T[k]/T_r[k]

!ALL IS MESSURE RESPECT TO T_r
T_d     = T_ref*T_r
!viscosity
!.....g/(cm*s)x10-6=Kg/(m*s)x10-7
mi_d    = dexp(A_mi*dlog(T_d) + B_mi/(T_d) + C_mi/(T_d)/(T_d) + D_mi) 

!Head conductinity
!.....W/(cm*K)x10-6
k_d     = dexp(A_k*dlog(T_d) + B_k/(T_d) + C_k/(T_d)/(T_d) + D_k ) 

!Head capacity 
! Cp is dimensionalees by the Ru, Cp=Cp/Ru
Cp_d    = A_cp + B_cp*T_d + C_cp*T_d**2.d0 + D_cp*T_d**3.d0 + E_cp*T_d**4.d0


s_s=0
e_s=N


!Dimesional variables
do j=s_s,e_s
        
        T_d     = T_ref*Tb(j)
        
        mib(j)  =  dexp(A_mi*dlog(T_d)  + B_mi/(T_d)  + C_mi/(T_d)/(T_d) + D_mi) 
        
        kb(j)   =  dexp( A_k*dlog(T_d)  + B_k /(T_d)  + C_k /(T_d)/(T_d) + D_k ) 
        
        Cpb(j)  =  (A_cp + B_cp*(T_d)   + C_cp*(T_d)**2.d0 + D_cp*(T_d)**3.d0 + E_cp*(T_d)**4.d0) 
!
END DO 

!Real value of cp, because is non dimesionalised by Ru
Cp_d  = Cp_d*Ru/Mw
Cpb   = Cpb *Ru/Mw
!Cp_d  = Cp_d*Ru
!Cpb   = Cpb *Ru
!write(*,*)'Cp_r',Cp_r*Ru/Mw1

END SUBROUTINE proper_1 

subroutine proper_2(mi_d,mib,k_d,kb,cp_d,cpb,Tb,T_ref,T_r,Mw,Ru,N)

implicit none
integer :: i,j
real(kind=ip)::T_b,T_ref,T_r,T_d
real(kind=ip)::A_mi,B_mi,C_mi,D_mi
real(kind=ip)::A_k , B_k ,C_k ,D_k 
real(kind=ip)::A_cp ,B_cp ,C_cp ,D_cp,E_cp
integer      ::N
real(kind=ip),dimension(0:N)::Tb
real(kind=ip),dimension(0:N)::mib,kb,cpb
real(kind=ip)::mi_d,k_d,cp_d
real(kind=ip)::Mw,Ru


CALL          visc_O2(A_mi,B_mi,C_mi,D_mi)
CALL          cond_O2(A_k ,B_k ,C_k ,D_k )
CALL  heatcapacity_O2(A_cp ,B_cp ,C_cp ,D_cp,E_cp,Mw)

!.....Dimensionaless equation for 
!.....temperature, using T_r T'=T[k]/T_r[k]

!ALL IS MESSURE RESPECT TO T_r
T_d     = T_ref*T_r
!viscosity
!.....g/(cm*s)x10-6=Kg/(m*s)x10-7
mi_d    = dexp(A_mi*dlog(T_d) + B_mi/(T_d) + C_mi/(T_d)/(T_d) + D_mi) 

!Head conductinity
!.....W/(cm*K)x10-6
k_d     = dexp(A_k*dlog(T_d) + B_k/(T_d) + C_k/(T_d)/(T_d) + D_k ) 

!Head capacity 
! Cp is dimensionalees by the Ru, Cp=Cp/Ru
Cp_d    = A_cp + B_cp*T_d + C_cp*T_d**2.d0 + D_cp*T_d**3.d0 + E_cp*T_d**4.d0


s_s=0
e_s=N


!Dimesional variables
do j=s_s,e_s

        T_d     = T_ref*Tb(j)
        
        mib(j)  =  dexp(A_mi*dlog(T_d)  + B_mi/(T_d)  + C_mi/(T_d)/(T_d) + D_mi) 
        
        kb(j)   =  dexp( A_k*dlog(T_d)  + B_k /(T_d)  + C_k /(T_d)/(T_d) + D_k ) 
        
        Cpb(j)  =  (A_cp + B_cp*(T_d)   + C_cp*(T_d)**2.d0 + D_cp*(T_d)**3.d0 + E_cp*(T_d)**4.d0) 

END DO 

!Real value of cp, because is non dimesionalised by Ru
Cp_d  = Cp_d*Ru/Mw
Cpb   = Cpb *Ru/Mw

!write(*,*)'Cp_r',Cp_r*Ru/Mw1

END SUBROUTINE proper_2 

subroutine proper_3(mi_d,mib,k_d,kb,cp_d,cpb,Tb,T_ref,T_r,Mw,Ru,N)

implicit none
integer :: i,j
real(kind=ip)::T_b,T_ref,T_r,T_d
real(kind=ip)::A_mi,B_mi,C_mi,D_mi
real(kind=ip)::A_k , B_k ,C_k ,D_k 
real(kind=ip)::A_cp ,B_cp ,C_cp ,D_cp,E_cp
integer      ::N
real(kind=ip),dimension(0:N)::Tb
real(kind=ip),dimension(0:N)::mib,kb,cpb
real(kind=ip)::mi_d,k_d,cp_d
real(kind=ip)::Mw,Ru


CALL          visc_N2(A_mi,B_mi,C_mi,D_mi)
CALL          cond_N2(A_k ,B_k ,C_k ,D_k )
CALL  heatcapacity_N2(A_cp ,B_cp ,C_cp ,D_cp,E_cp,Mw)

!.....Dimensionaless equation for 
!.....temperature, using T_r T'=T[k]/T_r[k]

!ALL IS MESSURE RESPECT TO T_r
T_d     = T_ref*T_r
!viscosity
!.....g/(cm*s)x10-6=Kg/(m*s)x10-7
mi_d    = dexp(A_mi*dlog(T_d) + B_mi/(T_d) + C_mi/(T_d)/(T_d) + D_mi) 

!Head conductinity
!.....W/(cm*K)x10-6
k_d     = dexp(A_k*dlog(T_d) + B_k/(T_d) + C_k/(T_d)/(T_d) + D_k ) 

!Head capacity 
! Cp is dimensionalees by the Ru, Cp=Cp/Ru
Cp_d    = A_cp + B_cp*T_d + C_cp*T_d**2.d0 + D_cp*T_d**3.d0 + E_cp*T_d**4.d0


s_s=0
e_s=N


!Dimesional variables
do j=s_s,e_s
        
        T_d     = T_ref*Tb(j)
        
        mib(j)  =  dexp(A_mi*dlog(T_d)  + B_mi/(T_d)  + C_mi/(T_d)/(T_d) + D_mi) 
        
        kb(j)   =  dexp( A_k*dlog(T_d)  + B_k /(T_d)  + C_k /(T_d)/(T_d) + D_k ) 
        
        Cpb(j)  =  (A_cp + B_cp*(T_d)   + C_cp*(T_d)**2.d0 + D_cp*(T_d)**3.d0 + E_cp*(T_d)**4.d0) 

END DO 

!Real value of cp, because is non dimesionalised by Ru
Cp_d  = Cp_d*Ru/Mw
Cpb   = Cpb *Ru/Mw

!write(*,*)'Cp_r',Cp_r*Ru/Mw1

END SUBROUTINE proper_3 

subroutine proper_4(mi_d,mib,k_d,kb,cp_d,cpb,Tb,T_ref,T_r,Mw,Ru,N)

implicit none
integer :: i,j
real(kind=ip)::T_b,T_ref,T_r,T_d
real(kind=ip)::A_mi,B_mi,C_mi,D_mi
real(kind=ip)::A_k , B_k ,C_k ,D_k 
real(kind=ip)::A_cp ,B_cp ,C_cp ,D_cp,E_cp
integer      ::N
real(kind=ip),dimension(0:N)::Tb
real(kind=ip),dimension(0:N)::mib,kb,cpb
real(kind=ip)::mi_d,k_d,cp_d
real(kind=ip)::Mw,Ru


CALL          visc_CH4(A_mi,B_mi,C_mi,D_mi)
CALL          cond_CH4(A_k ,B_k ,C_k ,D_k )
CALL  heatcapacity_CH4(A_cp ,B_cp ,C_cp ,D_cp,E_cp,Mw)

!.....Dimensionaless equation for 
!.....temperature, using T_r T'=T[k]/T_r[k]

!ALL IS MESSURE RESPECT TO T_r
T_d     = T_ref*T_r
!viscosity
!.....g/(cm*s)x10-6=Kg/(m*s)x10-7
mi_d    = dexp(A_mi*dlog(T_d) + B_mi/(T_d) + C_mi/(T_d)/(T_d) + D_mi) 

!Head conductinity
!.....W/(cm*K)x10-6
k_d     = dexp(A_k*dlog(T_d) + B_k/(T_d) + C_k/(T_d)/(T_d) + D_k ) 

!Head capacity 
! Cp is dimensionalees by the Ru, Cp=Cp/Ru
Cp_d    = A_cp + B_cp*T_d + C_cp*T_d**2.d0 + D_cp*T_d**3.d0 + E_cp*T_d**4.d0


s_s=0
e_s=N


!Dimesional variables
do j=s_s,e_s
        
        T_d     = T_ref*Tb(j)
        
        mib(j)  =  dexp(A_mi*dlog(T_d)  + B_mi/(T_d)  + C_mi/(T_d)/(T_d) + D_mi) 
        
        kb(j)   =  dexp( A_k*dlog(T_d)  + B_k /(T_d)  + C_k /(T_d)/(T_d) + D_k ) 
        
        Cpb(j)  =  (A_cp + B_cp*(T_d)   + C_cp*(T_d)**2.d0 + D_cp*(T_d)**3.d0 + E_cp*(T_d)**4.d0) 

END DO 

!Real value of cp, because is non dimesionalised by Ru
Cp_d  = Cp_d*Ru/Mw
Cpb   = Cpb *Ru/Mw

!write(*,*)'Cp_r',Cp_r*Ru/Mw1

END SUBROUTINE proper_4 

!###### NITROGEN_2 Ar[200-1000 K] ######
SUBROUTINE VISC_N2(A_d,B_d,C_d,D_d)

!IN-OUT 
!Coeficcient fitting polynome 
!paper:Transport Coefficients for the NASA Lewis Chemical Equilibrium Program
!author: Roger Svehla 1995
!A_MI
!B_MI
!C_MI
!D_MI
!Temperature range[200-1000][K]


IMPLICIT NONE
real(kind=ip)::A_d,B_d,C_d,D_d


	A_d	=	0.62526577D0 
	B_d	=	-31.779652D0 
	C_d	=	-1640.7983D0 
	D_d	=	 1.7454992D0 



END SUBROUTINE VISC_N2

SUBROUTINE Cond_N2(A_d,B_d,C_d,D_d)
!Temperature range[200-1000][K]
IMPLICIT NONE
real(kind=ip)::A_d,B_d,C_d,D_d

  
	A_d=	0.85372829D0 
	B_d=	 105.18665D0 
	C_d=	-12299.753D0 
	D_d=	0.48299104D0 

END SUBROUTINE cond_N2

SUBROUTINE heatcapacity_N2(A_d,B_d,C_d,D_d,E_d,Mw)
!Temperature range[1000-3000][K]
IMPLICIT NONE
real(kind=ip)::A_d,B_d,C_d,D_d,E_d,Mw

!200k<T<6000k

	A_d= 3.53100528D0 
	B_d=-1.23660987E-4 
        C_d=-5.02999437E-7 
	D_d= 2.43530612E-9 
	E_d=-1.40881235E-12 

        Mw = 28.01348D0
!T>1000k

!	A_d= 2.95257626D0 
!	B_d= 1.39690057E-3 
!        C_d=-4.92631691E-7 
!	D_d= 7.86010367E-11 
!	E_d=-4.60755321E-15 

!        Mw = 28.01348D0

END SUBROUTINE heatcapacity_N2

SUBROUTINE diffusion_N2(delta_df,epsilonk_df)

!IN-OUT 
!Coeficcient fitting polynome 
!paper:Viscosities and Thermal Conductivities 
!       of Gases at Hight Temperatures
!author: Roger Svehla 1962
!delta_df
!epsilonk_df

IMPLICIT NONE
real(kind=ip)::delta_df,epsilonk_df

        delta_df        = 3.798d0       
        epsilonk_df     =  71.4d0       

END SUBROUTINE diffusion_N2
!


SUBROUTINE visc_N(A_d,B_d,C_d,D_d)
!Temperature range[1000-3000][K]
IMPLICIT NONE
real(kind=ip)::A_d,B_d,C_d,D_d

	A_d=	0.82926975D0 
	B_d=	 405.82833D0 
	C_d=	-159002.42D0 
	D_d=	0.17740763D0 


END SUBROUTINE VISC_N

SUBROUTINE Cond_N(A_d,B_d,C_d,D_d)
!Temperature range[1000-3000][K]
IMPLICIT NONE
real(kind=ip)::A_d,B_d,C_d,D_d

	A_d=	0.82928303D0 
	B_d=	 405.77643D0 
	C_d=	-158950.37D0 
	D_d=	0.97751462D0 
END SUBROUTINE cond_N


SUBROUTINE heatcapacity_N(A_d,B_d,C_d,D_d,E_d,Mw)
!Temperature range[1000-3000][K]
IMPLICIT NONE
real(kind=ip)::A_d,B_d,C_d,D_d,E_d,Mw

!T<1000k

	A_d=2.5D0 
	B_d=0.D0 
	C_d=0.D0 
	D_d=0.D0 
	E_d=0.D0 

        Mw = 14.00674D0

END SUBROUTINE heatcapacity_N
!
SUBROUTINE diffusion_N(delta_df,epsilonk_df)

!IN-OUT 
!Coeficcient fitting polynome 
!paper:Viscosities and Thermal Conductivities 
!       of Gases at Hight Temperatures
!author: Roger Svehla 1962
!delta_df
!epsilonk_df

IMPLICIT NONE
real(kind=ip)::delta_df,epsilonk_df

        delta_df        = 3.798d0       
        epsilonk_df     =  71.4d0       

END SUBROUTINE diffusion_N
!
!
SUBROUTINE VISC_He(A_d,B_d,C_d,D_d)
!Temperature range[200-1000][K]
IMPLICIT NONE
real(kind=ip)::A_d,B_d,C_d,D_d
	A_d	=	0.75015944D0 
	B_d	=	 35.763243D0 
	C_d	=	-2212.1291D0 
	D_d	=	0.92126352D0 
END SUBROUTINE VISC_He

!###### HIDROGEN_2 Ar[200-1000 K] ######
SUBROUTINE VISC_H2(A_d,B_d,C_d,D_d)

!IN-OUT 
!Coeficcient fitting polynome 
!paper:Transport Coefficients for the NASA Lewis Chemical Equilibrium Program
!author: Roger Svehla 1995
!A_MI
!B_MI
!C_MI
!D_MI
!Temperature range[200-1000][K]


IMPLICIT NONE
real(kind=ip)::A_d,B_d,C_d,D_d


	A_d	=	0.74553182D0 
	B_d	=	43.555109D0 
	C_d	=      -3257.9340D0 
	D_d	=	0.13556243D0 

END SUBROUTINE VISC_H2

SUBROUTINE Cond_H2(A_d,B_d,C_d,D_d)
!Temperature range[200-1000][K]
IMPLICIT NONE
real(kind=ip)::A_d,B_d,C_d,D_d

  
	A_d=	1.0240124D0 
	B_d=	297.09752D0 
	C_d=   -31396.363D0 
	D_d=	1.0560824D0 

END SUBROUTINE cond_H2

SUBROUTINE heatcapacity_H2(A_d,B_d,C_d,D_d,E_d,Mw)
!Temperature range[1000-3000][K]
IMPLICIT NONE
real(kind=ip)::A_d,B_d,C_d,D_d,E_d,Mw

!200k<T<6000k

	A_d= 2.34433112D0
	B_d= 7.98052075E-3
        C_d=-1.94781510E-5
	D_d= 2.01572094E-8
	E_d=-7.37611761E-12

        Mw = 2.01588D0
!T>1000k

!	A_d= 2.95257626D0 
!	B_d= 1.39690057E-3 
!        C_d=-4.92631691E-7 
!	D_d= 7.86010367E-11 
!	E_d=-4.60755321E-15 

!        Mw = 28.01348D0

END SUBROUTINE heatcapacity_H2
!
SUBROUTINE diffusion_H2(delta_df,epsilonk_df)

!IN-OUT 
!Coeficcient fitting polynome 
!paper:Viscosities and Thermal Conductivities 
!       of Gases at Hight Temperatures
!author: Roger Svehla 1962
!delta_df
!epsilonk_df

IMPLICIT NONE
real(kind=ip)::delta_df,epsilonk_df

        delta_df        = 2.827d0       
        epsilonk_df     =  59.7d0       

END SUBROUTINE diffusion_H2
!
!###### OXIGEN_2 Ar[200-1000 K] ######
SUBROUTINE VISC_O2(A_d,B_d,C_d,D_d)

!IN-OUT 
!Coeficcient fitting polynome 
!paper:Transport Coefficients for the NASA Lewis Chemical Equilibrium Program
!author: Roger Svehla 1995
!A_MI
!B_MI
!C_MI
!D_MI
!Temperature range[200-1000][K]


IMPLICIT NONE
real(kind=ip)::A_d,B_d,C_d,D_d


	A_d  =    0.60916180D0 
	B_d  =   -52.244847D0 
	C_d  =   -599.74009D0 
	D_d  =    2.0410801D0 

END SUBROUTINE VISC_O2

SUBROUTINE COND_O2(A_d,B_d,C_d,D_d)
!Temperature range[200-1000][K]
IMPLICIT NONE
real(kind=ip)::A_d,B_d,C_d,D_d

  
	A_d  =    0.77238828D0 
	B_d  =    6.9293259D0 
	C_d  =   -5900.8518D0 
	D_d  =    1.2202965D0 

END SUBROUTINE COND_O2

SUBROUTINE heatcapacity_O2(A_d,B_d,C_d,D_d,E_d,Mw)
!Temperature range[200-6000][K]
IMPLICIT NONE
real(kind=ip)::A_d,B_d,C_d,D_d,E_d,Mw

!200k<T<6000k

	A_d=  3.78246636D0  
	B_d= -2.99673416E-03 
        C_d=  9.84730200E-06
	D_d= -9.68129608E-09  
	E_d=  3.24372836E-12

        Mw =  31.9988D0
!T>1000k

!	A_d= 2.95257626D0 
!	B_d= 1.39690057E-3 
!        C_d=-4.92631691E-7 
!	D_d= 7.86010367E-11 
!	E_d=-4.60755321E-15 

!        Mw = 28.01348D0

END SUBROUTINE heatcapacity_O2

SUBROUTINE diffusion_O2(delta_df,epsilonk_df)

!IN-OUT 
!Coeficcient fitting polynome 
!paper:Viscosities and Thermal Conductivities 
!       of Gases at Hight Temperatures
!author: Roger Svehla 1962
!delta_df
!epsilonk_df

IMPLICIT NONE
real(kind=ip)::delta_df,epsilonk_df

        delta_df        = 3.467d0       
        epsilonk_df     = 106.7d0       

END SUBROUTINE diffusion_O2

!###### OXIGEN_2 Ar[200-1000 K] ######
SUBROUTINE VISC_CH4(A_d,B_d,C_d,D_d)

!IN-OUT 
!Coeficcient fitting polynome 
!paper:Transport Coefficients for the NASA Lewis Chemical Equilibrium Program
!author: Roger Svehla 1995
!A_MI
!B_MI
!C_MI
!D_MI
!Temperature range[200-1000][K]


IMPLICIT NONE
real(kind=ip)::A_d,B_d,C_d,D_d


	A_d  =    0.57643622D0 
	B_d  =   -93.7040790D0 
	C_d  =     869.92395D0 
	D_d  =     1.7333347D0 

END SUBROUTINE VISC_CH4

SUBROUTINE COND_CH4(A_d,B_d,C_d,D_d)
!Temperature range[200-1000][K]
IMPLICIT NONE
real(kind=ip)::A_d,B_d,C_d,D_d

  
	A_d  =    0.77238828D0 
	B_d  =    6.9293259D0 
	C_d  =   -5900.8518D0 
	D_d  =    1.2202965D0 

END SUBROUTINE COND_CH4

SUBROUTINE heatcapacity_CH4(A_d,B_d,C_d,D_d,E_d,Mw)
!Temperature range[200-6000][K]
IMPLICIT NONE
real(kind=ip)::A_d,B_d,C_d,D_d,E_d,Mw

!200k<T<6000k

	A_d=  1.63552643D0  
	B_d=  1.00842795E-02 
        C_d= -3.36916254E-06
	D_d=  5.34958667E-10  
	E_d= -3.15518833E-14

        Mw =  16.04276D0
!T>1000k

!	A_d= 2.95257626D0 
!	B_d= 1.39690057E-3 
!        C_d=-4.92631691E-7 
!	D_d= 7.86010367E-11 
!	E_d=-4.60755321E-15 

!        Mw = 28.01348D0

END SUBROUTINE heatcapacity_CH4

SUBROUTINE diffusion_CH4(delta_df,epsilonk_df)

!IN-OUT 
!Coeficcient fitting polynome 
!paper:Viscosities and Thermal Conductivities 
!       of Gases at Hight Temperatures
!author: Roger Svehla 1962
!delta_df
!epsilonk_df

IMPLICIT NONE
real(kind=ip)::delta_df,epsilonk_df

        delta_df        = 3.758d0       
        epsilonk_df     = 148.6d0       

END SUBROUTINE diffusion_CH4

END MODULE

