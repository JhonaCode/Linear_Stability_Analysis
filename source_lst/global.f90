module global

      implicit none
      !Reference properties                            
      integer, parameter :: ip = selected_real_kind(15, 307)
      integer:: s_s,e_s 
      !size of vector
      !integer:: N_s 
      !Cp_0:head capacity
      !REAL(KIND=ip) :: l_0, a_0, rho_0!, T_0, mi_0,k_0,Cp_0      

      !Non dimensional numbers                           
      !REAL(KIND=IP) :: RE,PR,PE,ST,MA                    

      !-------------------------------------------------------------------
      ! Difussion Properties 
      !
      !!VISCOSITY VARIABLE PERTURBATION
      !real(kind=ip),allocatable,dimension(:,:) ::mi
      !!VISCOSITY VARIABLE BASE FLOW
      !real(kind=ip),allocatable,dimension(:,:) ::miba
      !!Conductivity VARIABLE
      !real(kind=ip),allocatable,dimension(:,:) ::ka
      !!Conductivity VARIABLE BASE FLOW
      !real(kind=ip),allocatable,dimension(:,:) ::kba
      !!Heat Capacity VARIABLE
      !real(kind=ip),allocatable,dimension(:,:) ::cp
      !!Heat Capacity VARIABLE BASE FLOW
      !real(kind=ip),allocatable,dimension(:,:) ::cpba
      

      !TEMPERATURE
      !TEMPERATURE
      !real(kind=ip),allocatable,dimension(:,:) ::Temp,Tempba
      !!TEMPERATURE Derivate
      !REAL(KIND=IP),ALLOCATABLE,DIMENSION(:,:) ::dTempdm,dTempdn               
      !REAL(KIND=IP),ALLOCATABLE,DIMENSION(:,:) ::d2Tempdm,d2Tempdn               

      !REAL(KIND=IP),ALLOCATABLE,DIMENSION(:,:) ::dTempbadm,dTempbadn               

      !REAL(KIND=IP),ALLOCATABLE,DIMENSION(:,:) ::d2Tempbadm,d2Tempbadn               

      !!Mc	= Mach center line velocity of the jet
      !!Bb 	= Thickness of the jet 
      !!Hb 	= Radio of the center line homegenuous velocity
      !TR 	= Reservatory temperature 
      !T0 	= Jet Center line temperature 

end module



