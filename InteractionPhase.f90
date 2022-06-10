 Program InteractionPhase  
!     Reading test trajectory
 implicit none 
	
	! First do loop that reads the file
	! define varables as real datatype (x, y, z) coordinates atom positions
	real x(50000), y(50000), z(50000) 
	! Defining variables integer type data
	integer i, natoms, molnum, eof, molidx, N, j, k, l 
	
	! i is?, natoms is number of atoms, molnum is number of molecules, eof is end of file? molidx is molecule index, N is?
	! defining the lab frame - at top for easy changing
	real Zlabx, Zlaby, Zlabz, Ylabx, Ylaby, Ylabz, Xlabx, Xlaby, Xlabz
	! Diagonal elements of the polarizability tensor of water at top for easy changing & safe proximity to pbc 
	real alpha(3), prox(3)
		
	character(len = 11) :: fname, mname, cryst, natom
	character(len = 11) :: molid, idx, cone
	
	! define varables as real datatype (x, y, z) periodic boundary condition lengths
	real pbcx, pbcy, pbcz, cutoff 
	! define varables as real datatype (x, y, z) periodic boundary condition angles
	real angx, angy, angz 
	real Pi, epsilon0
	
	! Second do loop set records position data, molecular frame, and centers of mass.
	! This loop set also calculates the molecular polarizabilities and PiM
	! (x,y,z) positions of O, H1, H2, and the dummy normal
	real Ox, Oy, Oz, H1Wx, H1Wy, H1Wz, H2Wx, H2Wy, H2Wz, dNormX, dNormY, dNormZ
	! H to O vectors; dummy then what we keep to file
	real O1HvecX, O1HvecY, O1HvecZ, O2HvecX, O2HvecY, O2HvecZ 
	real O1HX(50000), O1HY(50000), O1HZ(50000), O2HX(50000), O2HY(50000), O2HZ(50000)
	! cross product elements & normalized
	real NormX(50000), NormY(50000), NormZ(50000) 
	
	! Molecular Frame definitions 
	real zaxx, zaxy, zaxz, yaxx, yaxy, yaxz, xaxx, xaxy, xaxz !these are dummies
	real xaxisx(50000), xaxisy(50000), xaxisz(50000) !Molecular Frame x-axis
	real yaxisx(50000), yaxisy(50000), yaxisz(50000) !Molecular Frame y-axis
	real zaxisx(50000), zaxisy(50000), zaxisz(50000) !Molecular Frame z-axis
	
	! Matricies for Euler rotations and the final molecular polarizability tensor
	real, dimension(3,3) :: Euler, InvEuler, PolPrime, PiM
	! Euler angles for finding the Euler rotations
	real phi, theta, psi
	real, dimension(50000)	:: Polxx(50000), Polyy(50000),Polzz(50000) 
	real, dimension(50000)	:: Polxy(50000), Polxz(50000)
	real, dimension(50000)	:: Polyx(50000), Polyz(50000)
	real, dimension(50000)	:: Polzx(50000), Polzy(50000) ! PolPrime records
	! Center of Mass calculations for the dipole interaction 
	real CMass(3) !this is the dummy
	real, dimension(50000)	::	CMassx(50000), CMassy(50000), CMassz(50000) !saved positions
	common /com/  CMassx, CMassy, CMassz !saved positions attempt at getting them read in loop 3
	
	! Third do loop set checks pbc proximity and interaction polarizabilities and PiI
	! Dummy Variables and such for third do loop set
	real Holdmassx, Holdmassy, Holdmassz ! CoM position set for outer do loop
	real VaryMassx, VaryMassy, VaryMassz ! CoM position set for inner do loop
	real, dimension(3,3) 	:: HoldPol(3,3), VaryPol(3,3)! polarizability set for do loops
	real DistRadx, DistRady, DistRadz, DistRad
	real, dimension(3,3) 	:: dyad(3,3),  Trij(3,3), PiIntTen, PiI, PTensor
	logical Move, IntAct
	! for troubleshooting use
	real xout,yout,zout, radius, vecmag
	integer, pointer 	:: pt1, pt2
	integer, target 	:: t1, t2


!! Declaration section!!

! For the purposes of this code, the space fixed (lab) frame is defined as the the identity matrix.
	! x vector from I
	Xlabx = 1 
	Xlaby = 0 
	Xlabz = 0
	! y vector from I
	Ylabx = 0 
	Ylaby = 1 
	Ylabz = 0
	! z vector from I
	Zlabx = 0 
	Zlaby = 0
	Zlabz = 1  

! From J. Chem. Phys. 140, 054507 (2014); 
! https://doi.org/10.1063/1.4863696 #Actual numbers come from another paper, by Murphy in the 1970s
! This set can altered to study different substances from water.
	alpha(1)=1.528
	alpha(2)=1.415
	alpha(3)=1.468
	
	Pi = 3.1415962
	epsilon0 = 78.4
	cutoff = 3.5
	pbcx = 40
	pbcy = 40
	pbcz = 40
	
! For the purposes of this code, the proximity to the pbc is tested with prox(x,y,z) defining where we don't need to fuss
	prox(1) = pbcx - cutoff
	prox(2) = pbcy - cutoff
	prox(3) = pbcz - cutoff


! Opening files for read/write	
  open(unit=11, file='water_drop_test.pdb', status='old') !opening this pdb inputfile
  open(unit=12, file='IntPhase.out', status='new') !opens an output file to write to
  
! Setting the iteration to 0  
  i = 0
  PiM = 0 
  PiI = 0
  
! First do loop that reads in the data from the files.  
  do
	if (i .EQ. 0) then 
	 read(11,fmt=1000, iostat=eof) cryst, pbcx, pbcy, pbcz,angx, angy, angz 
1000     format(A6,3(3X,F6.3), 3(2X,F5.2))  
	 write(12,fmt=2000) cryst, pbcx, pbcy, pbcz, angx, angy, angz
2000     format(A6,3(3X,F6.3), 3(2X,F5.2))    
	else if (eof < 0 ) then !exits if failed
	  exit
!        else if (i > 5) then !is trying to make sure this is only a short amount for testing
!          exit 
	else 
	 read(11,fmt=3000, iostat=eof) natom, molidx, mname, molid, idx, molnum, x(i), y(i), z(i)
3000     format(A6, I5,1X, A4,1X, A3,1X ,A1, I4,4X, 3(F8.3))
!         print 3000, natom, molidx, mname, molid, idx, !writing it back out into the test.out file
!     & molnum, x(i), y(i), z(i)   
	 if (natom == 'CONECT') then
	 exit
	 end if 
	
	 write(12,FMT=4000) 'i = ', i, natom, molidx, mname, molid, idx, molnum, x(i), y(i), z(i)
 4000     format(A4, I6, 1X, A4, I7, A5, A4,1X,A1,1X, I4,4X, 3(F8.3))  
	 endif
	i = i + 1
  end do  
! End initial read/write of the atom positions  

! Begin second do loop 
! Setting the itterators  
  N = i - 1 
  write(*,*) 'N = ', N 
  j = 1
  
  do i=5,N !actual start of the pdb coord is at position i=5
	if (mod(i-1,3) .eq. 0) then  
	! Writes the coordinates to screen 
		! write(*,*) 'Molecule ', i 
		! write(*,*) 'Coordinates:' 
		! write(*,*)  x(i), y(i), z(i)
		
	! Provides the labels for each atom in context so we can run calcs
		Ox=x(i)
		Oy=y(i)
		Oz=z(i)
		H1Wx=x(i-2)
		H1Wy=y(i-2)
		H1Wz=z(i-2)
		H2Wx=x(i-1)
		H2Wy=y(i-1)
		H2Wz=z(i-1)
		
!	Calculates a dummy atom normal to the plane defined by the H-O-H bonds and those vectors		
	call Dummy_Atom(j, Ox, Oy, Oz, H1Wx, H1Wy, H1Wz, H2Wx, H2Wy, H2Wz, dNormX, dNormY, dNormZ, &
		O1HvecX, O1HvecY, O1HvecZ, O2HvecX, O2HvecY, O2HvecZ)

!	H1W to O vector		
		O1HX(j) = O1HvecX
		O1HY(j) = O1HvecY
		O1HZ(j) = O1HvecZ
	! write(12,FMT=5000) 'H1W to O vector', j, O1HX(j), O1HY(j), O1HZ(j)
! 5000		format(A15,1X, I4,4X, 3(F8.3))  

!	H2W to O vector
		O2HX(j) = O2HvecX 
		O2HY(j) = O2HvecY
		O2HZ(j) = O2HvecZ
	! write(12,FMT=6000) 'H2W to O vector', j, O2HX(j), O2HY(j), O2HZ(j)
! 6000	format(A15,1X, I4,4X, 3(F8.3))  

!	Normal vector to the molecular plane		
		NormX(j)=dNormX 
		NormY(j)=dNormY 
		NormZ(j)=dNormZ
	! write(12,FMT=7000) 'Dummy Normal', j, NormX(j), NormY(j), NormZ(j)
! 7000	format(A12,4X, I4,4X, 3(F8.3))  
		
!	Molecular Frame part
!	This section records the specifcs of the molecular frame for the Euler rotation.
	call Molframe(j, H1Wx, H1Wy, H1Wz, H2Wx, H2Wy, H2Wz,  &
	 	O1HvecX, O1HvecY, O1HvecZ, O2HvecX, O2HvecY, O2HvecZ, &
		zaxx, zaxy, zaxz, xaxx, xaxy, xaxz)
		
		xaxisx(j) = xaxx
		xaxisy(j) = xaxy
		xaxisz(j) = xaxz
	! write(12,FMT=8000) 'Molframe X-axis', j, xaxisx(j), xaxisy(j), xaxisz(j)
! 8000	format(A15,1X, I4,4X, 3(F8.3))

		yaxisx(j) = NormX(j)
		yaxisy(j) = NormY(j)
		yaxisz(j) = NormZ(j)
	! write(12,FMT=9000) 'Molframe Y-axis', j, yaxisx(j), yaxisy(j), yaxisz(j)
! 9000     format(A15,1X, I4,4X, 3(F8.3))  


		zaxisx(j) = zaxx 
		zaxisy(j) = zaxy
		zaxisz(j) = zaxz
	! write(12,FMT=10000) 'Molframe Z-axis', j, zaxisx(j), zaxisy(j), zaxisz(j)
! 10000	format(A15,1X, I4,4X, 3(F8.3))  
		
		
!	My own coding of the Euler angles and rotation
	! We save this part in the output file
	call EulerAng(xaxx, xaxy, xaxz, Zlabx, Zlaby, Zlabz, zaxx, zaxy, zaxz, &
					Xlabx, Xlaby, Xlabz, phi, theta, psi)
	! write(12,FMT=11000) 'Euler Angles', j, phi, theta, psi
! 11000		format(A15,1X, I4,4X, 3(F8.3))  					
	
	
	! My own coding of the Euler Rotation matrix and its inverse.
	! This bit gets discarded for output file as unnecessary to keep.
	call EulerRot(phi, theta, psi, Euler, InvEuler)
	
	! My own coding of PolarRot = InvEuler*polar*Euler. 
	! This is the bit we save in the output file because it sums later.
	call PolarRot (Euler, InvEuler, alpha, PolPrime)
		Polxx(j) = PolPrime(1,1)
		Polxy(j) = PolPrime(1,2)
		Polxz(j) = PolPrime(1,3)
		Polyx(j) = PolPrime(2,1)
		Polyy(j) = PolPrime(2,2)
		Polyz(j) = PolPrime(2,3)
		Polzx(j) = PolPrime(3,1)
		Polzy(j) = PolPrime(3,2)
		Polzz(j) = PolPrime(3,3)
	write(*,*) 'PolPrime' , PolPrime
	write(12,FMT=12000) 'PolPrime', j, PolPrime
12000	format(A17,1X, I4,12X, 3(F15.3),/, 34x,3(F15.3),/, 34x,3(F15.3),/)  	
	
	call CenterMass(H1Wx, H1Wy, H1Wz, H2Wx, H2Wy, H2Wz, Ox, Oy, Oz, CMass)
		CMassx(j) = CMass(1)
		CMassy(j) = CMass(2)
		CMassz(j) = CMass(3)
	! write(*,*) 'Center of Mass', CMass
	! write(12,FMT=13000) 'Center of Mass for', j, CMassx(j), CMassy(j), CMassz(j)
! 13000	format(A15,1X, I4,4X, 3(F8.3)) 

!	This gives us the sum of the molecular polarizabilities	Pi(M)
	PiM = PiM + PolPrime


!	Next, to find the sum of the interaction polarizabilities Pi(I)
	j = j + 1	
		
	endif
end do

	
! Third set of do loops to check the interactions 
! Setting the iteration to 0  
! nullify(pt1)
! nullify(pt2)
k = 1	

	do k = 1, j-1 ! start at one and go through each molecule
	
		! Hold is for the i in the summation, Vary is for the j such that i/=j
		! Center of Masses called in for dummy variables
		Holdmassx= CMassx(k)
		Holdmassy= CMassy(k)
		Holdmassz= CMassz(k)
		! write(*,*) 'Center of mass held array', k, CMassx(k), CMassy(k), CMassz(k)
		! write(*,*) 'Center of mass held read in', Holdmassx, Holdmassy, Holdmassz
		
		! ! Polarization  element for each CoM(k)
		HoldPol(1,1) = Polxx(k)
		HoldPol(1,2) = Polxy(k) 
		HoldPol(1,3) = Polxz(k)
		
		HoldPol(2,1) = Polyx(k)
		HoldPol(2,2) = Polyy(k)
		HoldPol(2,3) = Polyz(k)
		
		HoldPol(3,1) = Polzx(k)
		HoldPol(3,2) = Polzy(k)
		HoldPol(3,3) = Polzz(k)
		
		! There should be few molecules next to the pbc, thus run as normal.
		!call PBCcheck(prox, Holdmassx, Holdmassy, Holdmassz, Move)
			! if (Move .eqv. .true.) then
			! write(*,*) 'prox', k, Move
			!else
						
			do l = k+1, j ! all post-k molecules
						
				! ! Dummy variables Center of Masses that will vary wrt k
				VaryMassx= CMassx(l)
				VaryMassy= CMassy(l)
				VaryMassz= CMassz(l)
				! !write(*,*) 'Center of mass varried array', k, l, CMassx(l), CMassy(l), CMassz(l)
				! !write(*,*) 'Center of mass varried read in', VaryMassx,VaryMassy,VaryMassz
								
					
				call IntDistCheck (cutoff, Holdmassx,Holdmassy,Holdmassz, &
							VaryMassx,VaryMassy,VaryMassz, IntAct, DistRadx,DistRady,DistRadz, DistRad)
					if (IntAct .eqv. .true.) then
					! write(*,*) 'Interacting?', IntAct
					
					call IntTensor (Pi, epsilon0, DistRadx,DistRady,DistRadz,DistRad, Trij)
						write(*,*) 'Interaction i to j', k, l
					! write(12,FMT=16000) 'Interaction Tensor', k, l, Trij
! 16000					format(A17,1X, I4,4X,I4,4X, 3(F15.3),/, 34x,3(F15.3),/, 34x,3(F15.3),/) 
					
						! !!Polarization element for each CoM(k)				
						VaryPol(1,1) = Polxx(l)
						VaryPol(1,2) = Polxy(l) 
						VaryPol(1,3) = Polxz(l)
						
						VaryPol(2,1) = Polyx(l)
						VaryPol(2,2) = Polyy(l)
						VaryPol(2,3) = Polyz(l)
						
						VaryPol(3,1) = Polzx(l)
						VaryPol(3,2) = Polzy(l)
						VaryPol(3,3) = Polzz(l)
						
						call PiInt (HoldPol, VaryPol, Trij, PiIntTen)
							! write(*,*) 'Polarization Pi I', PiIntTen
						write(12,FMT=17000) 'PiIntensor',k, l, PiIntTen
	17000					format(A17,1X, I4,4X,I4,4X, 3(F15.3),/, 34x,3(F15.3),/, 34x,3(F15.3),/)		
						
						PiI = PiI + PiIntTen
						write(*,*) 'Pi I is now:', PiI
					
					else
						cycle !skips back should interaction /=  true
					endif
				
							
			end do
			! else !this part will happen once i have fully troubleshot it
				 ! cycle
			! endif
			
		end do
	
!	Eventually, I'll put this at the end of file so that Pi(M) is next to Pi(I)	
	write(12,FMT=14000) 'PiM, Molecular', PiM
14000		format(A17,17X, 3(F15.3),/, 34x,3(F15.3),/, 34x,3(F15.3),/)
	
	write(12,FMT=18000) 'PiI, Interaction',  PiI
18000		format(A17, 17X, 3(F15.3),/, 34x,3(F15.3),/, 34x,3(F15.3),/)

PTensor = PiM +PiI

write(12,FMT=19000) 'Collective Pi', PTensor
19000		format(A17, 17X, 3(F15.3),/, 34x,3(F15.3),/, 34x,3(F15.3))
	
end program InteractionPhase

!Subroutines
Subroutine IntDistCheck (cutoff, Holdmassx,Holdmassy,Holdmassz, &
							VaryMassx,VaryMassy,VaryMassz, IntAct, DistRadx,DistRady,DistRadz, DistRad)
	! input arguments
	implicit none
	real				:: Holdmassx, Holdmassy, Holdmassz
	real				:: VaryMassx, VaryMassy, VaryMassz
	! dummies
	real 				:: DistRadx,DistRady,DistRadz
	real				:: cutoff, radius, DistRad
	logical				:: IntAct
	
	call DipoleDist ( Holdmassx,Holdmassy,Holdmassz, &
			Varymassx, Varymassy, Varymassz, DistRadx,DistRady,DistRadz, DistRad)	
			
	radius = DistRad
	!write (*,*) 'radius', radius,'cutoff', cutoff
	IntAct = (radius .le. cutoff) 
	

end Subroutine IntDistCheck
Subroutine PBCcheck (prox, Holdmassx, Holdmassy, Holdmassz, Move)
	! input arguments
	implicit none
	real, dimension(3)	:: prox 
	real				:: Holdmassx, Holdmassy, Holdmassz
	logical				:: Ansx, Ansy, Ansz, Move
	
	Ansx = (prox(1) .le. Holdmassx) 
	
	Ansy = (prox(2) .le. Holdmassy) 
	
	Ansz = (prox(3) .le. Holdmassz) 

	Move =((Ansx .eqv. .true.) .or. (Ansy .eqv. .true.) .or. (Ansz .eqv. .true.)) 

end Subroutine PBCcheck
Subroutine PiInt (HoldPol, VaryPol, Trij, PiIntTen)
	! input arguments
	implicit none
	real, dimension(3,3) :: HoldPol, VaryPol, Trij
	real, dimension(3,3) :: PiIntTen, temp
	
	temp = matmul(Trij,VaryPol)
	PiIntTen = matmul(HoldPol, temp)!*Trij*VaryPol

end Subroutine PiInt		
Subroutine IntTensor (Pi, epsilon0, DistRadx,DistRady,DistRadz, DistRad, Trij)
	! input arguments
	implicit none      
	real, intent(in)	 :: DistRadx,DistRady,DistRadz
	real, intent(in)	 :: Pi, epsilon0,  DistRad
	! dummy arguments    
	real				 :: mag, permiability, scalar
	real 				 :: VecMag
	real, dimension(3,3) :: Unittensor(3,3), dyad(3,3), temp(3,3), ScaMatrix(3,3)
	real, dimension(3,3) :: temp2(3,3), Trij(3,3)
	
	scalar = 3
	
	Unittensor(1,1) = 1
	Unittensor(2,2) = 1
	Unittensor(3,3) = 1
	Unittensor(1,2) = 0
	Unittensor(1,3) = 0
	Unittensor(2,1) = 0
	Unittensor(2,3) = 0
	Unittensor(3,1) = 0
	Unittensor(3,2) = 0
	
	call Dyadic(DistRadx,DistRady,DistRadz,DistRadx,DistRady,DistRadz,dyad)
	
	
	call ScalarMult(scalar, dyad, ScaMatrix) 
	temp = ScaMatrix
	temp2 = temp - Unittensor
	
	mag = DistRad**3
	
	permiability = 1/(4*pi*epsilon0*mag)
	
	call ScalarMult(permiability, temp2, ScaMatrix)
	Trij = ScaMatrix

end Subroutine
Subroutine Dyadic(ax1,ay2,az3,bx1,by2,bz3,dyad)
	! Dyadic tensor product aka outer product circle cross symbol 
	! input arguments
	implicit none      
	real, intent(in)	 :: ax1,ay2,az3
	real, intent(in)	 :: bx1,by2,bz3
	! dummy arguments        
	real, dimension(3,3) :: dyad(3,3)
	
	dyad(1,1)= ax1*bx1
	dyad(1,2)= ax1*by2
	dyad(1,3)= ax1*bz3
	
	dyad(2,1)= ay2*bx1
	dyad(2,2)= ay2*by2
	dyad(2,3)= ay2*bz3
	
	dyad(3,1)= az3*bx1
	dyad(3,2)= az3*by2
	dyad(3,3)= az3*bz3
	
end Subroutine Dyadic
Subroutine DipoleDist ( Holdmassx,Holdmassy,Holdmassz, &
			Varymassx, Varymassy, Varymassz, DistRadx,DistRady,DistRadz, DistRad)	
	! input arguments
	implicit none
	real	:: Holdmassx,Holdmassy,Holdmassz
	real	:: Varymassx, Varymassy, Varymassz
	! dummy arguments
	real	:: DistRadx, DistRady, DistRadz, DistRad
	real 	:: VecMag, dradout
	
	call dist (Holdmassx,Holdmassy,Holdmassz, Varymassx, Varymassy, Varymassz, &
			DistRadx,DistRady,DistRadz)
				
	DistRad = VecMag(DistRadx, DistRady, DistRadz)
		
		
end Subroutine DipoleDist
Subroutine CenterMass (x1,y1,z1, x2,y2,z2, x3,y3,z3, CMass) 
	! input arguments
	implicit none
	real, intent(in)	:: x1,y1,z1, x2,y2,z2, x3,y3,z3
	! dummy arguments
	real, dimension(3)	:: CMass
	real				:: Hmass, Omass
	
	Omass = 15.999		! https://www.webelements.com/oxygen/ just a fave online table
	Hmass = 1.008 		! https://www.webelements.com/hydrogen/
	
	CMass(1) = (1/(2* Hmass + Omass)) * (Hmass*(x1+x2)+ Omass*x3)
	CMass(2) = (1/(2* Hmass + Omass)) * (Hmass*(y1+y2)+ Omass*y3)
	CMass(3) = (1/(2* Hmass + Omass)) * (Hmass*(z1+z2)+ Omass*z3)
	
end subroutine CenterMass
Subroutine PolarRot(Euler, InvEuler, alpha, PolPrime)
	! input arguments
	implicit none      
	real, intent(in) 	 :: Euler(3,3), InvEuler(3,3)
	! dummy arguments        
	real, dimension(3)   :: alpha(3)
	real, dimension(3,3) :: temp(3,3), Polar(3,3), PolPrime(3,3)
	
	! Define the polarizability tensor of the molecule.
	Polar(1,1)=alpha(1)
	Polar(2,2)=alpha(2)
	Polar(3,3)=alpha(3)
	Polar(1,2)=0
	Polar(1,3)=0
	Polar(2,1)=0
	Polar(2,3)=0
	Polar(3,1)=0
	Polar(3,2)=0

	temp = matmul(polar,Euler)
	PolPrime = matmul(InvEuler, temp)!*(polar*Euler)
	
end subroutine PolarRot
Subroutine EulerRot(phi, theta, psi, Euler, InvEuler)

	! input arguments
	implicit none      
	real, intent(in) :: phi, theta, psi
	! dummy arguments        
	real, dimension(3,3) :: Euler(3,3), InvEuler(3,3)
	
	Euler(1,1)= cos(phi) * cos(psi) - sin(phi) * cos(theta) * sin(psi)
	Euler(1,2)= -cos(phi) * sin(psi) - sin(phi) * cos(theta) * cos(psi)
	Euler(1,3)= sin(phi) * sin(theta)
	Euler(2,1)= sin(phi)* cos(psi) + cos(phi) * cos(theta) * sin(psi)
	Euler(2,2)= -sin(phi)* sin(psi) + cos(phi) * cos(theta) * cos(psi)
	Euler(2,3)= -cos(phi) * sin(theta)
	Euler(3,1)= sin(theta) * sin(psi)
	Euler(3,2)= sin(theta) * cos(psi)
	Euler(3,3)= cos(theta)
	
	!InvEuler is also Euler(transpose) by definition and orthonomality.
	InvEuler(1,1)= cos(phi) * cos(psi) - sin(phi) * cos(theta) * sin(psi)
	InvEuler(2,1)= -cos(phi) * sin(psi) - sin(phi) * cos(theta) * cos(psi)
	InvEuler(3,1)= sin(phi) * sin(theta)
	InvEuler(1,2)= sin(phi)* cos(psi) + cos(phi) * cos(theta) * sin(psi)
	InvEuler(2,2)= -sin(phi)* sin(psi) + cos(phi) * cos(theta) * cos(psi)
	InvEuler(3,2)= -cos(phi) * sin(theta)
	InvEuler(1,2)= sin(theta) * sin(psi)
	InvEuler(2,3)= sin(theta) * cos(psi)
	InvEuler(3,3)= cos(theta)
	
end Subroutine EulerRot	
Subroutine EulerAng (xaxx, xaxy, xaxz, Zlabx, Zlaby, Zlabz, zaxx, zaxy, zaxz, &
					Xlabx, Xlaby, Xlabz, phi, theta, psi)
	! input arguments
	implicit none      
	real, intent(in)  :: Zlabx, Zlaby, Zlabz, zaxx, zaxy, zaxz
	real, intent(in)  :: Xlabx, Xlaby, Xlabz, xaxx, xaxy, xaxz
	! dummy arguments        
	real			  :: angle, xout, yout, zout, nodex, nodey, nodez
	real, intent(out) :: phi, theta, psi
	
	! first find the line of nodes, this is the vector n = Zlab x Zmolframe
	call Crossproduct(Zlabx, Zlaby, Zlabz, zaxx, zaxy, zaxz, xout, yout, zout)
	nodex = xout
	nodey = yout
	nodez = zout
	
	! phi is the angle between the Xmolframe and the vector n.
	phi = angle (xaxx, xaxy, xaxz, nodex, nodey, nodez)
		
	! theta is the angle between the space (lab) frame Zlab and the Zmolframe 
	theta = angle (Zlabx, Zlabz, Zlabz, zaxx, zaxy, zaxz)
	if (theta .eq. 0) then
	  theta = 0.001
	end if !theta cannot equal zero or we lose linearity
	
	! psi is the angle between the line of nodes, n, and Xmolframe.
	psi = angle (nodex, nodey, nodez, Xlabx, Xlaby, Xlabz)
	
end subroutine EulerAng	
Subroutine Molframe(j, H1Wx, H1Wy, H1Wz, H2Wx, H2Wy, H2Wz, &
		O1HvecX, O1HvecY, O1HvecZ, O2HvecX, O2HvecY, O2HvecZ, &
		zaxx, zaxy, zaxz, xaxx, xaxy, xaxz)

	! The oxygen atom is set as the origin for the molecular frame.
	! x-axis
	! This is calculated through finding the line that bisects the triangle
	! formed by the H-O-H bond through the oxygen in the molecular frame. 
	! need O1HVec and O2Hvec and the VecMag
	call xaxis(O1HvecX, O1HvecY, O1HvecZ, O2HvecX, O2HvecY, O2HvecZ, xout2, yout2, zout2)
	xaxx = xout2
	xaxy = yout2
	xaxz = zout2

	! y-axis
	! The y-axis is orthogonal to the z-axis and x-axis such that the Oxygen
	! lone pairs in the classic vespr model of water are in the XY plain 
	! with the C2v symmetry preserved.
		! yaxx = dNormX
		! yaxy = dNormY 
		! yaxz = dNormZ
	
	
	! z-axis 
	! This is calculated by finding the vector between the hydrogens and
	! translating the z variable to the oxygen position in the molecular frame. 
	call zaxis(H1Wx, H1Wy, H1Wz, H2Wx, H2Wy, H2Wz, xaxz, xout,yout,zout)
	zaxx = xout
	zaxy = yout
	zaxz = zout

end Subroutine Molframe
Subroutine Dummy_Atom(j, Ox, Oy, Oz, H1Wx, H1Wy, H1Wz, H2Wx, H2Wy, H2Wz, dNormX, dNormY, dNormZ, &
				O1HvecX, O1HvecY, O1HvecZ, O2HvecX, O2HvecY, O2HvecZ)
	
	integer, intent(in) :: j
	! Declaration for internal variables
	! Variables for the creation of the dummy atom on the normal vector
	real O1HvecX, O1HvecY, O1HvecZ !dipole vector/bond lengths along the H1W to O bond
	real O2HvecX, O2HvecY, O2HvecZ !dipole vector/bond lengths along the H2W to O bond
	real magO1H, magO2H, NormFact  !magnitudes of the two bond lengths & normalization factor
	real CrosX, CrosY, CrosZ !cross product elements & normalized

!	Distance from O to H1W by element - function call to follow
	call dist (Ox, Oy, Oz, H1Wx, H1Wy, H1Wz, O1HvecX, O1HvecY, O1HvecZ)
	
!	Distance from O to H2W by element - function call to follow	
	call dist (Ox, Oy, Oz, H2Wx, H2Wy, H2Wz, O2HvecX, O2HvecY, O2HvecZ )
		
!	Magnitude of the vector from O to H1W - function call to follow 		
	magO1H = VecMag(O1HvecX, O1HvecY, O1HvecZ)
!	Magnitude of the vector from O to H2W - function call to follow
	magO2H = VecMag(O2HvecX, O2HvecY, O2HvecZ)
!	Normalization factor using (normal)sin (pi/2) = 1 = ((O to H1W)x(O to H2W))/(magO1H * magO2H)
	NormFact = magO1H * magO2H	
	
!	Cross product of O-H1W and O-H2W by element - function call to follow will include the normalization
	call Crossproduct(O1HvecX, O1HvecY, O1HvecZ, O2HvecX, O2HvecY, O2HvecZ, xout, yout, zout)
	CrosX = xout
	CrosY = yout
	CrosZ = zout
	
!	Normalization of the crossproduct so that our normal dummy atom is 1 angstrom up
	dNormX = CrosX/NormFact
	dNormY = CrosY/NormFact
	dNormZ = CrosZ/NormFact
	
end subroutine Dummy_Atom
Subroutine zaxis(x1,y1,z1,x2,y2,z2,xaxz, xout,yout,zout)

	! input arguments
	implicit none      
	real, intent(in) :: x1,y1,z1,x2,y2,z2, xaxz
	! dummy arguments        
	real :: magHH, HHx, HHy, HHz, VecMag
	real :: zdirx, zdiry, zdirz
	real ::	xout,yout,zout

	call dist(x1,y1,z1,x2,y2,z2, HHx, HHy, HHz)
	zdirx = HHx
	zdiry = HHy
	zdirz = HHz 
	magHH = VecMag(zdirx, zdiry, zdirz)

	xout = HHx / magHH
	yout = HHy / magHH
	zout = (HHz / magHH) + xaxz !translation to Ow
		
end Subroutine zaxis
Subroutine xaxis(x1,y1,z1,x2,y2,z2, xout, yout, zout)

	! input arguments
	implicit none      
	real, intent(in) :: x1,y1,z1,x2,y2,z2
	! dummy arguments        
	real :: HHtheta, normO1H, normO2H
	real :: divisor, vecmag, angle
	real ::	xout, yout, zout

	HHtheta = angle(x1,y1,z1,x2,y2,z2)
	normO1H = VecMag(x1,y1,z1)
	normO2H = VecMag(x2,y2,z2)
	divisor = (2 * (normO1H * normO2H)) / HHtheta

	xout = x1 / divisor
	yout = y1 / divisor
	zout = z1 / divisor
			
end Subroutine xaxis
Subroutine dist(x1,y1,z1,x2,y2,z2,xout,yout,zout)

	! input arguments
	implicit none      
	real, intent(in) :: x1,y1,z1,x2,y2,z2
	! dummy arguments        
	real ::  xout,yout,zout
	
	xout = x1-x2
	yout = y1-y2
	zout = z1-z2
	
end Subroutine dist	
Subroutine ScalarMult (scalar, Matrix,ScaMatrix)
	! input arguments
	implicit none      
	real, intent(in)	 :: scalar
	real, intent(in)	 :: Matrix(3,3)
	! dummy arguments        
	real, dimension(3,3) :: ScaMatrix(3,3)
	
	ScaMatrix(1,1)= scalar*Matrix(1,1)
	ScaMatrix(1,2)= scalar*Matrix(1,2)
	ScaMatrix(1,3)= scalar*Matrix(1,3)
	
	ScaMatrix(2,1)= scalar*Matrix(2,1)
	ScaMatrix(2,2)= scalar*Matrix(2,2)
	ScaMatrix(2,3)= scalar*Matrix(2,3)
	
	ScaMatrix(3,1)= scalar*Matrix(3,1)
	ScaMatrix(3,2)= scalar*Matrix(3,2)
	ScaMatrix(3,3)= scalar*Matrix(3,3)
	
end Subroutine ScalarMult
Subroutine Crossproduct(x1,y1,z1, x2,y2,z2, xout, yout, zout)
	! input arguments
	implicit none      
	real, intent(in)  :: x1,y1,z1,x2,y2,z2
	! dummy arguments
	!real			  :: CrosX, CrosY, CrosZ
	real			  :: CrossProdX, CrossProdY, CrossProdZ
	real, intent(out) :: xout, yout, zout

	xout = CrossProdX(y1, z1, y2, z2)	
	yout = CrossProdY(x1, z1, x2, z2)
	zout = CrossProdZ(x1, y1, x2, y2)
 
 end Subroutine Crossproduct

 !! Functions 
 ! real function SquareMult(MatrixA(3,3), MatrixB(3,3)) result (MatrixC(3,3))
	! ! input arguments
	! implicit none      
	! real, intent(in)	 :: scalar
	! real, intent(in)	 :: MatrixA(3,3), MatrixB(3,3)
	! ! dummy arguments        
	! real, dimension(3,3) :: MatrixC(3,3)
	
	! MatrixC(1,1)= MatrixA(1,1) * MatrixB(1,1)
	! MatrixC(1,2)= MatrixA(1,2) * MatrixB(3,3)
	! MatrixC(1,3)= MatrixA(1,3) * MatrixB(3,3)
	
	! MatrixC(2,1)= MatrixA(2,1) * MatrixB(1,1)
	! MatrixC(2,2)= MatrixA(2,2) * MatrixB(3,3)
	! MatrixC(2,3)= MatrixA(2,3) * MatrixB(3,3)
	
	! MatrixC(3,1)= MatrixA(3,1) * MatrixB(1,1)
	! MatrixC(3,2)= MatrixA(3,2) * MatrixB(3,3)
	! MatrixC(3,3)= MatrixA(3,3) * MatrixB(3,3)
	
 ! end function SquareMult
real function CrossProdX(y1, z1, y2, z2)  result (CrosX)

	! input arguments     
	implicit none      
	real, intent(in)  :: y1, z1, y2, z2
	! dummy arguments        
	!real, intent(out) :: CrosX
	 
	CrosX = y1 * z2 - z1 * y2
	!	CrosY = z1 * x2 - x1 * z2
	!	CrosZ = x1 * y2 - y1 * x2

end function CrossProdX	
real function CrossProdY(x1, z1, x2, z2)  result (CrosY)

	! input arguments     
	implicit none      
	real, intent(in)  :: x1, z1, x2, z2
	! dummy arguments        
	!real, intent(out) :: CrosY    
	 
	!	CrosX = y1 * z2 - z1 * y2
	CrosY = z1 * x2 - x1 * z2
	!	CrosZ = x1 * y2 - y1 * x2

end function CrossProdY	
real function CrossProdZ(x1, y1, x2, y2)  result (CrosZ)

	! input arguments     
	implicit none      
	real, intent(in)  :: x1, y1, x2, y2
	! dummy arguments        
	!real, intent(out) ::  CrosZ   
	 
	!	CrosX = y1 * z2 - z1 * y2
	!	CrosY = z1 * x2 - x1 * z2
	CrosZ = x1 * y2 - y1 * x2

end function CrossProdZ	
real function angle(x1,y1,z1,x2,y2,z2) result (theta)

	! function result     
	implicit none      
	real, intent(in) :: x1,y1,z1,x2,y2,z2
	! dummy arguments    
	! real :: theta
	real ::  mag1, mag2, mag, doot 
	
	mag1 = sqrt (x1**2 + y1**2 + z1**2)
	mag2 = sqrt (x2**2 + y2**2 + z2**2)
	mag = mag1 * mag2 
	doot = (x1*x2+y1*y2+z1*z2)
	!acos is inverse cosine in radians!
	theta = acos(doot/mag)
	
end function angle
real function VecMag(x1,y1,z1)  result (MagVec)

	! function result     
	implicit none      
	real :: x1,y1,z1
	! dummy arguments         
 
	MagVec = sqrt (x1**2 + y1**2 + z1**2)
	
end function VecMag	
real function dotprod(x1,y1,z1,x2,y2,z2) result(dotscalar)
	! function result     
	implicit none      
	real, intent(in) :: x1,y1,z1,x2,y2,z2
	! dummy arguments        
	! real ::  dotscalar

	dotscalar = (x1*x2+y1*y2+z1*z2)
	
end function dotprod	

