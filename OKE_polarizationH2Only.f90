 Program OKE_polarizationH2Only
!   Program calculates the collective polarizability of a water box at one instant in time.
!	
!	Still to be done: 
!		1. read in of pbc from user input on screen / slurm 
!		2. automate the count of atoms and molecules
!			
 implicit none 
	! as good practice, specify the level of precision for the numbers.
	! First do loop that reads the file
	! define varables as real datatype (x, y, z) coordinates atom positions
	real, dimension(:), allocatable ::  x, y, z 
	! Defining variables integer type data
	integer i, natoms, molnum, eof, molidx
	integer N, j, k, l 
	integer t, timesteps, atom, totatom
	
	! i is?, natoms is number of atoms, molnum is number of molecules, eof is end of file? molidx is molecule index, N is?
	! defining the lab frame - at top for easy changing
	real Zlabx, Zlaby, Zlabz, Ylabx, Ylaby, Ylabz, Xlabx, Xlaby, Xlabz
	! Diagonal elements of the polarizability tensor of water at top for easy changing & safe proximity to pbc 
	real alpha(3), prox(3)
	real Gammasqr, TCFdenom
		
	character(len = 11) :: fname, mname, cryst, natom
	character(len = 11) :: molid, idx, cone
	
	! define varables as real datatype (x, y, z) periodic boundary condition lengths
	real pbcx, pbcy, pbcz, cutoff 
	! define varables as real(kind=dp)dp datatype (x, y, z) periodic boundary condition angles
	real angx, angy, angz 
	real Pi, epsilon0
	
	! Second do loop set records position data, molecular frame, and centers of mass.
	! This loop set also calculates the molecular polarizabilities and PiM
	! (x,y,z) positions of O, H1, H2, and the dummy normal
	real Ox, Oy, Oz, H1Wx, H1Wy, H1Wz, H2Wx, H2Wy, H2Wz, dNormX, dNormY, dNormZ
	! H to O vectors; dummy then what we keep to file
	real O1HvecX, O1HvecY, O1HvecZ, O2HvecX, O2HvecY, O2HvecZ 
	real O1HX, O1HY, O1HZ, O2HX, O2HY, O2HZ
	! cross product elements & normalized
	real NormX, NormY, NormZ 
	
	! Molecular Frame definitions 
	real zaxx, zaxy, zaxz, yaxx, yaxy, yaxz, xaxx, xaxy, xaxz !these are dummies
	real xaxisx, xaxisy, xaxisz !Molecular Frame x-axis
	real yaxisx, yaxisy, yaxisz !Molecular Frame y-axis
	real zaxisx, zaxisy, zaxisz !Molecular Frame z-axis
	
	! Matricies for Euler rotations and the final molecular polarizability tensor
	real, dimension(3,3):: Euler, InvEuler, PolPrime, PiM
	!Check the number of digits used with this definition of "real precision" numbers
	
	! Euler angles for finding the Euler rotations
	real phi, theta, psi
	
	real, dimension(:), allocatable :: Polxx, Polxy, Polxz
	real, dimension(:), allocatable :: Polyx, Polyy, Polyz
	real, dimension(:), allocatable :: Polzx, Polzy, Polzz ! PolPrime records
	! Center of Mass calculations for the dipole interaction 
	real CMass(3) !this is the dummy
	real, dimension(:), allocatable :: 	CMassx, CMassy, CMassz !saved positions
	!common /com/  CMassx, CMassy, CMassz !saved positions attempt at getting them read in loop 3
	
	! Third do loop set checks pbc proximity and interaction polarizabilities and PiI
	! Dummy Variables and such for third do loop set
	real HoldMassX, HoldMassY, HoldMassZ ! CoM position set for outer do loop
	real VaryMassX, VaryMassy, VaryMassZ ! CoM position set for inner do loop
	real NVarx, NVary, NVarz !New VaryMass should there be potential interactions across the pbc
	real, dimension(3,3) 	:: HoldPol(3,3), VaryPol(3,3)! polarizability set for do loops
	real DistRadx, DistRady, DistRadz, DistRad
	real, dimension(3,3) 	:: dyad(3,3),  Trij(3,3), PiIntTen, PiI, PTensor
	logical Move, IntAct
	! for troubleshooting use
	real xout,yout,zout, radius, vecmag
	


!! Declaration section!!

! For the initial hardcode using Nitesh's pdb file, these are true. They will be replaced with a dynamic subroutine later.
	timesteps = 2002-1 !count starts at zero, not 1
	molnum = 900 !totatom atoms, use this for array allocation
	totatom = molnum*3 !total number of atoms is 3 (H-O-H) times 900 molecules
	
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
! A second version of this should be run to observe the differences of the Murphy values vs the Huiszoon 1986 values
	alpha(1)=1.528
	alpha(2)=1.415
	alpha(3)=1.468
	
	Gammasqr = 0.5 * ( (alpha(1)-alpha(2))**2 +(alpha(1)-alpha(3))**2 +(alpha(2)-alpha(3))**2 )
	TCFdenom = (Gammasqr * molnum) / 15
	
!	Pi = 3.1415962
!	epsilon0 = 78.4 !not needed for this model
	cutoff = 3.5
!	pbcx = 40 !hardcode not needed for this model
!	pbcy = 40 !hardcode not needed for this model
!	pbcz = 40 !hardcode not needed for this model
		
! Opening files for read/write	
  open(unit=11, file='water_box_3A_nvt.pbc', status='old') !opening this pdb inputfile
  open(unit=12, file='OKE_polarizationH2Only.out', status='new') !opens an output file to write to
    open(unit=13, file='Pi_xz_at_T.out', status='new') !opens an output file to write to individual Pi-xz_at_T values wrt timestep
	!open(unit=14, file='PolarizationTCF.out', status='new') !opens an output file to write to TCF
    !open(unit=15, file='PiM_AllocTestsReadIn.out', status='new')
  
! Setting the iteration to 0 and initializing the Pis as zero
  i = 0
  j = 0
  t = 0
  PiM(:,:) = 0.d0 !real number and need to do the digits
  PiI(:,:) = 0.d0 !real number and need to do the digits
  PTensor(:,:) = 0.d0 !real number and need to do the digits
write(12,FMT=14000) 'PiM, Clear1', PiM !this was a check to be sure PiM was clear at start
write(12,FMT=14000) 'PiI, Clear1', PiI !this was a check to be sure PiI was clear at start
write(12,FMT=14000) 'PTensor, Clear1', PTensor !this was a check to be sure PTensor was clear at start

allocate (x(totatom))
allocate (y(totatom))
allocate (z(totatom))
allocate (Polxx(totatom))
allocate (Polxy(totatom))
allocate (Polxz(totatom))
allocate (Polyx(totatom))
allocate (Polyy(totatom))
allocate (Polyz(totatom))
allocate (Polzx(totatom))
allocate (Polzy(totatom))
allocate (Polzz(totatom))
allocate (CMassx(totatom))
allocate (CMassy(totatom))
allocate (CMassz(totatom))

read(11,fmt=1000, iostat=eof) cryst, pbcx, pbcy, pbcz,angx, angy, angz 
	1000     format(A6,3(3X,F6.3), 3(2X,F5.2))  
		  ! write(12,fmt=2000) cryst, pbcx, pbcy, pbcz, angx, angy, angz
	! 2000     format(A6,3(3X,F6.3), 3(2X,F5.2))    
	
 do t = 0, timesteps, 1 !do loop itterating over snapshots wrt time
		
	do i = 0, (totatom*(1+timesteps) +1)! Read in do loop that reads in the data from the files.  
	
		if (i .gt. 0) then 
		 
		read(11,fmt=3000, iostat=eof) natom, molidx, mname, molid, idx, molnum, x(i), y(i), z(i)
	3000    format(A6, I5,1X, A4,1X, A3,1X ,A1, I4,4X, 3(F8.3))
			if (eof < 0 ) then !exits if failed
				exit
			else if (natom .EQ. 'END') then !is trying to make sure this is only a short amount for testing
				exit
			else
		! print 3000, natom, molidx, mname, molid, idx, & !writing it back out into the test.out file
			! molnum, x(i), y(i), z(i)   
				
		! write(15,FMT=4000) 'i = ', i, natom, molidx, mname, molid, idx, molnum, x(i), y(i), z(i)
	! 4000     format(A4, I6, 1X, A4, I7, A5, A4,1X,A1,1X, I4,4X, 3(F8.3))  
								
		end if
		
	 
	 endif
	end do  ! end do Read in do loop that reads in the data from the files.  

	! For the purposes of this code, the proximity to the pbc is tested with prox(x,y,z) defining where we don't need to fuss
		prox(1) = pbcx - cutoff
		prox(2) = pbcy - cutoff
		prox(3) = pbcz - cutoff

	!need to start a new do loop to set the snapshot as timestep = [0, t]

	! Begin second do loop to find Pi(M)
	! Setting the itterators  
	
	
	N = i  ! i is itterating over atoms and is line position
	!write(*,*) 'N = ', N 
	! atom = 0

	
	  ! in total, j will equal molnum
	  
!do j = 1 , 900, 1 !not needed for this work

	j = 1
	do atom = 4, N, 3 !actual start of the pdb coord is at position i=2
		 ! j itterates over molecules, 3 atoms (i) per H-O-H molecule 
		
		! if (j .eq. 901) then
			! exit
		! else 
		!if (mod(atom+2,3) .eq. 0) then  
		! Writes the coordinates to screen 
			  ! write(*,*) 'Molecule ', atom 
			  ! write(*,*) 'Coordinates:' 
			  ! write(*,*)  x(atom), y(atom), z(atom)
		! write(15,FMT=5000) 'if read in check ', j, atom, x(atom), y(atom), z(atom)
! 5000		format(A15,1X, I4, 1X, I4,4X, 3(F8.3))  

		! Provides the labels for each atom in context so we can run calcs
			Ox=x(atom-2)
			Oy=y(atom-2)
			Oz=z(atom-2)
			H1Wx=x(atom-1)
			H1Wy=y(atom-1)
			H1Wz=z(atom-1)
			H2Wx=x(atom)
			H2Wy=y(atom)
			H2Wz=z(atom)
		
! Calculates a dummy atom normal to the plane defined by the H-O-H bonds and those vectors		
		call Dummy_Atom(j, Ox, Oy, Oz, H1Wx, H1Wy, H1Wz, H2Wx, H2Wy, H2Wz, dNormX, dNormY, dNormZ, &
			O1HvecX, O1HvecY, O1HvecZ, O2HvecX, O2HvecY, O2HvecZ)

			! H1W to O vector		
				O1HX = O1HvecX
				O1HY = O1HvecY
				O1HZ = O1HvecZ
			 ! write(12,FMT=5000) 'H1W to O vector', j, O1HX, O1HY, O1HZ
		! 5000		format(A15,1X, I4,4X, 3(F8.3))  

		
			! H2W to O vector
				O2HX = O2HvecX 
				O2HY = O2HvecY
				O2HZ = O2HvecZ
			! write(12,FMT=6000) 'H2W to O vector', j, O2HX, O2HY, O2HZ
		! 6000	format(A15,1X, I4,4X, 3(F8.3))  

			! Normal vector to the molecular plane		
				NormX = dNormX 
				NormY = dNormY 
				NormZ = dNormZ
			! write(12,FMT=7000) 'Dummy Normal', j, NormX, NormY, NormZ
		! 7000	format(A12,4X, I4,4X, 3(F8.3))  
			
! Molecular Frame part
		! This section records the specifcs of the molecular frame for the Euler rotation.
		call Molframe(j, H1Wx, H1Wy, H1Wz, H2Wx, H2Wy, H2Wz,  &
			O1HvecX, O1HvecY, O1HvecZ, O2HvecX, O2HvecY, O2HvecZ, &
			zaxx, zaxy, zaxz, xaxx, xaxy, xaxz)
			
				xaxisx = xaxx
				xaxisy = xaxy
				xaxisz = xaxz
			! write(12,FMT=8000) 'Molframe X-axis', j, xaxisx, xaxisy, xaxisz
		! 8000	format(A15,1X, I4,4X, 3(F8.3))

				yaxisx = NormX
				yaxisy = NormY
				yaxisz = NormZ
			! write(12,FMT=9000) 'Molframe Y-axis', j, yaxisx, yaxisy, yaxisz
		! 9000     format(A15,1X, I4,4X, 3(F8.3))  

				zaxisx = zaxx 
				zaxisy = zaxy
				zaxisz = zaxz
			! write(12,FMT=10000) 'Molframe Z-axis', j, zaxisx, zaxisy, zaxisz
		! 10000	format(A15,1X, I4,4X, 3(F8.3))  
			
! This subroutine calculates the Euler angles from the molecular frame
		call EulerAng(xaxx, xaxy, xaxz, Zlabx, Zlaby, Zlabz, zaxx, zaxy, zaxz, &
						Xlabx, Xlaby, Xlabz, phi, theta, psi)
			! write(15,FMT=11000) 'Euler Angles', t, j, phi, theta, psi
		! 11000		format(A15,1X, I4, 1X, I4,4X, 3(F8.3))  					
		
	! Logic check for invalid values (i.e. NaN errors), when present these are tossed and not used. 
			if (phi/phi .ne. 1) then
				! PolPrime = 0
				! j = j + 1
				exit
			elseif  (theta/theta .ne. 1) then
				! PolPrime = 0
				! j = j + 1
				exit
			elseif (psi/psi .ne. 1) then
				! PolPrime = 0
				! j = j + 1
				exit
			endif
		
	! This subroutine creates the Euler rotation matrix and its inverse	
		call EulerRot(phi, theta, psi, Euler, InvEuler)


! This subroutine actually calculates our individual molecular polarizabilities, Pi(M) 	
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
		! write(*,*) 'PolPrime' , PolPrime
			! write(15,FMT=12000) 'PolPrime', j, PolPrime
		! 12000	format(A17,1X, I4,12X, 3(F15.3),/, 34x,3(F15.3),/, 34x,3(F15.3),/)  
		
! Now the individual PiM is done, so alpha(i) in the lit		
		
! This subroutine just calculates the molecular centers of mass
		call CenterMass(H1Wx, H1Wy, H1Wz, H2Wx, H2Wy, H2Wz, Ox, Oy, Oz, CMass)
			CMassx(j) = CMass(1)
			CMassy(j) = CMass(2)
			CMassz(j) = CMass(3)
			! ! write(*,*) 'Center of Mass', CMass
			! write(15,FMT=13000) 'Center of Mass for', j, CMassx(j), CMassy(j), CMassz(j)
		 13000	format(A15,1X, I4,4X, 3(F8.3)) 

! This gives us the running sum of the molecular polarizabilities	Pi(M)
		PiM = PiM + PolPrime
			

	j = j + 1
	end do !atom = 2+t*(2701), N	
!end do !j = 1 , 900, 1 !not needed for this work	
	
! Interaction Polarization, Pi(I)	
! Third set of do loops to check the interactions 
	PiI = 0
	!	write(12,FMT=18000) 'PiI, Int clear 1',  PiI
	! 18000		format(A17, 17X, 3(F15.3),/, 34x,3(F15.3),/, 34x,3(F15.3),/)

	! Calcluates the total collective polarization tensor at this timestep
	PTensor = 0

	!write(12,FMT=19000) 'Col Pi cleared', PTensor
	19000		format(A17, 17X, 3(F15.5),/, 34x,3(F15.5),/, 34x,3(F15.5),/)
	k = 1	

		do k = 1, j-1 ! start at one and go through each molecule, i in the lit.
		
	! Hold is for the i in the summation, Vary is for the j such that i/=j in the lit
			! Center of Masses called in for dummy variables
			HoldMassX= CMassx(k)
			HoldMassY= CMassy(k)
			HoldMassZ= CMassz(k)
			 ! write(15,FMT=13000) 'Center of mass held array', k, CMassx(k), CMassy(k), CMassz(k)
			 ! write(15,FMT=13000) 'Center of mass held read in',k, HoldMassX, HoldMassY, HoldMassZ
			
		! Polarization  element for each CoM(k)
			HoldPol(1,1) = Polxx(k)
			HoldPol(1,2) = Polxy(k) 
			HoldPol(1,3) = Polxz(k)
			
			HoldPol(2,1) = Polyx(k)
			HoldPol(2,2) = Polyy(k)
			HoldPol(2,3) = Polyz(k)
			
			HoldPol(3,1) = Polzx(k)
			HoldPol(3,2) = Polzy(k)
			HoldPol(3,3) = Polzz(k)
			
	! There should be few molecules next to the pbc, the check for those is in the next if statement.
			call PBCcheck1(prox, HoldMassX, HoldMassY, HoldMassZ, Move)
										
				do l = k+1, j ! all post-k molecules, l is j in the lit for our Trij tensor
							
					! Dummy variables Center of Masses that will vary wrt k
					VaryMassX= CMassx(l)
					VaryMassy= CMassy(l)
					VaryMassZ= CMassz(l)
					! write(15,FMT=5000) 'Center of mass varried array',k, l, CMassx(l), CMassy(l), CMassz(l)
					! write(15,FMT=5000) 'Center of mass varried read in',k, l, VaryMassX,VaryMassy,VaryMassZ
					
	! Checking for potential interactions across the pbc
		! Only translate the molecules we must.	
					if (Move .eqv. .true.) then 
						call PBCcheck2 (prox,pbcx, pbcy, pbcz, HoldMassX, HoldMassY, HoldMassZ, &
									VaryMassX, VaryMassy, VaryMassZ, NVarx, NVary, NVarz)
							VaryMassX= NVarx
							VaryMassy= NVary
							VaryMassZ= NVarz
					endif		
					
	! Checks if k and l are close enough to interact					
					call IntDistCheck (cutoff, HoldMassX,HoldMassY,HoldMassZ, &
								VaryMassX,VaryMassy,VaryMassZ, IntAct, DistRadx,DistRady,DistRadz, DistRad)
						if (IntAct .eqv. .true.) then 
						! write(*,*) 'Interacting?', IntAct
						
						! Calculates Trij, the interaction tensor
						call IntTensor (Pi, epsilon0, DistRadx,DistRady,DistRadz,DistRad, Trij)
						! write(*,*) 'Interaction i to j', k, l
						! write(12,FMT=16000) 'Interaction Tensor', k, l, Trij
	! 16000					format(A17,1X, I4,4X,I4,4X, 3(F15.3),/, 34x,3(F15.3),/, 34x,3(F15.3),/) 
						
		! Polarization element for each CoM(k) that interacts			
							VaryPol(1,1) = Polxx(l)
							VaryPol(1,2) = Polxy(l) 
							VaryPol(1,3) = Polxz(l)
							
							VaryPol(2,1) = Polyx(l)
							VaryPol(2,2) = Polyy(l)
							VaryPol(2,3) = Polyz(l)
							
							VaryPol(3,1) = Polzx(l)
							VaryPol(3,2) = Polzy(l)
							VaryPol(3,3) = Polzz(l)
							
		! calls the subroutine to calculate the individual k to l polarization
							call PiInt (HoldPol, VaryPol, Trij, PiIntTen)
								! write(*,*) 'Polarization Pi I', PiIntTen
								! write(12,FMT=17000) 'PiIntensor',k, l, PiIntTen
! 17000								format(A17,1X, I4,4X,I4,4X, 3(F15.3),/, 34x,3(F15.3),/, 34x,3(F15.3),/)		
								
	! running summation of the net interactions between k and l
							PiI = PiI + PiIntTen
							! write(*,*) 'Pi I is now:', PiI
						
						else
							cycle !skips back should interaction /=  true
						endif
												
				end do ! finishing inner loop over l
							
			! end do ! finishing inner loop over k !not needed here
		
		
	j = j + 1
	
	end do ! j <= molnum, (j .lt. 900) 
	

	write(12, FMT=20000) 'timestep', t
20000		format(A17, 5X, I10.1,/)

	write(12,FMT=14000) 'PiM, Molecular', PiM
14000		format(A17,17X, 3(F17.3),/, 34x,3(F17.3),/, 34x,3(F17.3),/)
		
	write(12,FMT=18000) 'PiI, Interaction',  PiI
18000		format(A17, 17X, 3(F15.3),/, 34x,3(F15.3),/, 34x,3(F15.3),/)

! Calcluates the total collective polarization tensor at this timestep
PTensor = PiM +PiI

	write(12,FMT=19000) 'Collective Pi', PTensor
! 19000		format(A17, 17X, 3(F15.5),/, 34x,3(F15.5),/, 34x,3(F15.5),/)

	write(13,FMT=21000) 'Pi tensor element Pi-xz at timestep', PTensor(3,1), t
21000		format(A35, 10X, (F17.5), I10)
	
	! write(*,*) 'timestep', t
	! write(*,FMT=19000) 'Collective Pi', PTensor
  
	
 
  
!t = t + 1
!x(:)=0
!y(:)=0
!z(:)=0
!deallocate (x, y, z)	
!write(*,*) 'x, y, z, t are', x, y, z, t	
PiM(:,:) = 0.d0 !real number and need to do the digits
PiI(:,:) = 0.d0 !real number and need to do the digits
!write(12,FMT=14000) 'PiM, M-Clear', PiM !this was a check to be sure PiM was clearing between timesteps

end do !(ends the do while t <=timesteps)
!	Closing all files that were opened before the timesteps as no longer being read/written
	close(11)
	close(12)
	close(13)

	! open(unit=13, file='Pi_xz_at_T.out', status='old') !opens an output file to write to individual Pi-xz_at_T values wrt timestep
	! open(unit=14, file='Pi_xz_TCF.out', status='new')

	! read(13,FMT=21000) 'Pi tensor element Pi-xz at timestep', PTensor(3,1), t
!21000		format(A35, 10X, (F17.5), I10)
	! do p = 1, timesteps, 1

		! call TCF (Pten1, Pten2, TCFdenom, PtatT)
		! write(13,FMT=21000) 'Psi tensor element Psi-xz at timestep', PTensor(3,1), t
! 21000		format(A35, 10X, (F17.5), I10)
	! end do

! close(13)
! close(14)

end program OKE_polarizationH2Only

!Subroutines
Subroutine TCF
	! input arguments
	implicit none
	real				:: Pten1, Pten2, TCFdenom
	real				:: PtatT
	
	PtatT = (Pten1 * Pten2) / TCFdenom

end Subroutine TCF
Subroutine IntDistCheck (cutoff, HoldMassX,HoldMassY,HoldMassZ, &
							VaryMassX,VaryMassy,VaryMassZ, IntAct, DistRadx,DistRady,DistRadz, DistRad)
	! input arguments
	implicit none
	real				:: HoldMassX, HoldMassY, HoldMassZ
	real				:: VaryMassX, VaryMassy, VaryMassZ
	! dummies
	real 				:: DistRadx,DistRady,DistRadz
	real				:: cutoff, radius, DistRad
	logical				:: IntAct
	
	call DipoleDist ( HoldMassX,HoldMassY,HoldMassZ, &
			VaryMassX, VaryMassY, VaryMassZ, DistRadx,DistRady,DistRadz, DistRad)	
			
	radius = DistRad
	!write (*,*) 'radius', radius,'cutoff', cutoff
	IntAct = (radius .le. cutoff) 
	

end Subroutine IntDistCheck
Subroutine PBCcheck2 (prox,pbcx, pbcy, pbcz, HoldMassX, HoldMassY, HoldMassZ, &
							VaryMassX, VaryMassy, VaryMassZ, NVarx, NVary, NVarz)
	! input arguments
	implicit none
	real, dimension(3)	:: prox 
	real				:: pbcx, pbcy, pbcz
	real				:: HoldMassX, HoldMassY, HoldMassZ
	real				:: VaryMassX, VaryMassy, VaryMassZ 
	real				:: NVarx, NVary, NVarz 
	
	if (VaryMassX .lt. abs(prox(1) - HoldMassX)) then
		NVarx = VaryMassX + pbcx !translates on x-axis
		else
		NVarx = VaryMassX
	endif 
	
	if (VaryMassy .lt. abs(prox(2) - HoldMassY)) then
		NVary = VaryMassy + pbcy !translates on y-axis
		else
		NVary = VaryMassy
	endif 
	
	if (VaryMassZ .lt. abs(prox(3) - HoldMassZ)) then
		NVarz = VaryMassZ + pbcz !translates on z-axis 
		else
		NVarz = VaryMassZ
	endif 
	
end Subroutine PBCcheck2	
Subroutine PBCcheck1 (prox, HoldMassX, HoldMassY, HoldMassZ, Move)
	! input arguments
	implicit none
	real, dimension(3)	:: prox 
	real				:: HoldMassX, HoldMassY, HoldMassZ
	logical				:: Ansx, Ansy, Ansz, Move
	
	Ansx = (prox(1) .le. HoldMassX) 
	
	Ansy = (prox(2) .le. HoldMassY) 
	
	Ansz = (prox(3) .le. HoldMassZ) 

	Move =((Ansx .eqv. .true.) .or. (Ansy .eqv. .true.) .or. (Ansz .eqv. .true.)) 

end Subroutine PBCcheck1
Subroutine PiInt (HoldPol, VaryPol, Trij, PiIntTen)
	! input arguments
	implicit none
	real, dimension(3,3) :: HoldPol, VaryPol, Trij
	real, dimension(3,3) :: PiIntTen, temp
	
	temp = matmul(Trij,VaryPol)
	PiIntTen = matmul(HoldPol, temp)

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
	
	permiability = 1/mag !(4*pi*epsilon0*mag)
	
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
Subroutine DipoleDist ( HoldMassX,HoldMassY,HoldMassZ, &
			VaryMassX, VaryMassY, VaryMassZ, DistRadx,DistRady,DistRadz, DistRad)	
	! input arguments
	implicit none
	real	:: HoldMassX,HoldMassY,HoldMassZ
	real	:: VaryMassX, VaryMassY, VaryMassZ
	! dummy arguments
	real	:: DistRadx, DistRady, DistRadz, DistRad
	real 	:: VecMag, dradout
	
	call dist (HoldMassX,HoldMassY,HoldMassZ, VaryMassX, VaryMassY, VaryMassZ, &
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

