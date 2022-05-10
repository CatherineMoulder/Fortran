 Program EulerRotCoor  
!     Reading test trajectory
  implicit none 
	real x(50000), y(50000), z(50000) !define varables as real datatype (x, y, z) coordinates
	integer i, natoms, molnum, eof, molidx, N, j !Defining variables integer type data
	!i is?, natoms is number of atoms, molnum is number of molecules, eof is end of file? molidx is molecule index, N is?
  
	real Zlabx, Zlaby, Zlabz, Ylabx, Ylaby, Ylabz, Xlabx, Xlaby, Xlabz
	character(len = 11) :: fname, mname, cryst, natom
	character(len = 11) :: molid, idx, cone
	real pbex, pbey, pbez !define varables as real datatype (x, y, z) Euler coordinates in the labframe
	real angx, angy, angz !define varables as real datatype Euler angles in the labframe
	real Ox, Oy, Oz, H1Wx, H1Wy, H1Wz, H2Wx, H2Wy, H2Wz, dNormX, dNormY, dNormZ
	real O1HvecX, O1HvecY, O1HvecZ, O2HvecX, O2HvecY, O2HvecZ !H to O vectors
	real O1HX(50000), O1HY(50000), O1HZ(50000), O2HX(50000), O2HY(50000), O2HZ(50000)
	real NormX(50000), NormY(50000), NormZ(50000) !cross product elements & normalized
	real zaxx, zaxy, zaxz, yaxx, yaxy, yaxz, xaxx, xaxy, xaxz !these are dummies
	real xaxisx(50000), xaxisy(50000), xaxisz(50000) !Molecular Frame x-axis
	real yaxisx(50000), yaxisy(50000), yaxisz(50000) !Molecular Frame y-axis
	real zaxisx(50000), zaxisy(50000), zaxisz(50000) !Molecular Frame z-axis
	real alpha(3), CMass(3)!     Diagonal elements of the polarizability tensor of water 
	real, dimension(3,3) :: Euler, InvEuler, PolPrime, Ptensor
	real phi, theta, psi
	

!From J. Chem. Phys. 140, 054507 (2014); 
!https://doi.org/10.1063/1.4863696 #Actual numbers come from another paper, by Murphy in the 1970s
!This set can altered to study different substances from water.
	alpha(1)=1.528
	alpha(2)=1.415
	alpha(3)=1.468

!For the purposes of this code, the space fixed (lab) frame is defined as the the identity matrix.
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
	
  open(unit=11, file='water_drop_test.pdb', status='old') !opening this pdb inputfile
  open(unit=12, file='Cmass.out', status='new') !opens an output file to write to
  i = 0
  do
	if (i .EQ. 0) then 
	 read(11,fmt=1000, iostat=eof) cryst, pbex, pbey, pbez,angx, angy, angz 
1000     format(A6,3(3X,F6.3), 3(2X,F5.2))  
	 write(12,fmt=2000) cryst, pbex, pbey, pbez, angx, angy, angz
2000     format(A6,3(3X,F6.3), 3(2X,F5.2))    
	else if (eof < 0 ) then !exits if failed
	  exit
!        else if (i > 5) then 
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
  
  N = i - 1 
  write(*,*) 'N = ', N 
  j=1
  do i=5,N !actual start of the pdb coord is at position i=5
	if (mod(i-1,3) .eq. 0) then  
		write(*,*) 'Molecule ', i 
		write(*,*) 'Coordinates:' 
		write(*,*)  x(i), y(i), z(i)
	
		Ox=x(i)
		Oy=y(i)
		Oz=z(i)
		H1Wx=x(i-2)
		H1Wy=y(i-2)
		H1Wz=z(i-2)
		H2Wx=x(i-1)
		H2Wy=y(i-1)
		H2Wz=z(i-1)
		
		call Dummy_Atom(j, Ox, Oy, Oz, H1Wx, H1Wy, H1Wz, H2Wx, H2Wy, H2Wz, dNormX, dNormY, dNormZ, &
		O1HvecX, O1HvecY, O1HvecZ, O2HvecX, O2HvecY, O2HvecZ)

!	H1W to O vector		
		O1HX(j) = O1HvecX
		O1HY(j) = O1HvecY
		O1HZ(j) = O1HvecZ
		write(12,FMT=5000) 'H1W to O vector', j, O1HX(j), O1HY(j), O1HZ(j)
5000		format(A15,1X, I4,4X, 3(F8.3))  

!	H2W to O vector
		O2HX(j) = O2HvecX 
		O2HY(j) = O2HvecY
		O2HZ(j) = O2HvecZ
		write(12,FMT=6000) 'H2W to O vector', j, O2HX(j), O2HY(j), O2HZ(j)
6000		format(A15,1X, I4,4X, 3(F8.3))  

!	Normal vector to the molecular plane		
		NormX(j)=dNormX 
		NormY(j)=dNormY 
		NormZ(j)=dNormZ
		write(12,FMT=7000) 'Dummy Normal', j, NormX(j), NormY(j), NormZ(j)
7000		format(A12,4X, I4,4X, 3(F8.3))  
		
!	Molecular Frame part
!	This section records the specifcs of the molecular frame for the Euler rotation.
		call Molframe(j, H1Wx, H1Wy, H1Wz, H2Wx, H2Wy, H2Wz,  &
	 	O1HvecX, O1HvecY, O1HvecZ, O2HvecX, O2HvecY, O2HvecZ, &
		zaxx, zaxy, zaxz, xaxx, xaxy, xaxz)
		
		xaxisx(j) = xaxx
		xaxisy(j) = xaxy
		xaxisz(j) = xaxz
		write(12,FMT=8000) 'Molframe X-axis', j, xaxx, xaxy, xaxz
8000		format(A15,1X, I4,4X, 3(F8.3))

		yaxisx(j) = NormX(j)
		yaxisy(j) = NormY(j)
		yaxisz(j) = NormZ(j)
		write(12,FMT=9000) 'Molframe Y-axis', j, yaxisx(j), yaxisy(j), yaxisz(j)
9000     format(A15,1X, I4,4X, 3(F8.3))  


		zaxisx(j) = zaxx 
		zaxisy(j) = zaxy
		zaxisz(j) = zaxz
		write(12,FMT=10000) 'Molframe Z-axis', j, zaxx, zaxy, zaxz
10000		format(A15,1X, I4,4X, 3(F8.3))  
		
		
!	My own coding of the Euler angles 
	! We save this part in the output file
	call EulerAng(xaxx, xaxy, xaxz, Zlabx, Zlaby, Zlabz, zaxx, zaxy, zaxz, &
					Xlabx, Xlaby, Xlabz, phi, theta, psi)
	write(12,FMT=11000) 'Euler Angles', j, phi, theta, psi
11000		format(A15,1X, I4,4X, 3(F8.3))  					
	
	
	! My own coding of the Euler Rotation matrix and its inverse.
	! This bit gets discarded for output file as unnecessary to keep.
	call EulerRot(phi, theta, psi, Euler, InvEuler)
	
	! My own coding of PolarRot = InvEuler*polar*Euler. 
	! This is the bit we save in the output file because it sums later.
	call PolarRot (Euler, InvEuler, alpha, PolPrime)
	write(12,FMT=12000) 'PolPrime Diagonal', j, PolPrime(1,1), PolPrime(2,2), PolPrime(3,3)
12000		format(A17,1X, I4,4X, 3(F8.3))  	
	
	call CenterMass(H1Wx, H1Wy, H1Wz, H2Wx, H2Wy, H2Wz, Ox, Oy, Oz, CMass)
	write(12,FMT=13000) 'Center of Mass for', j, CMass
13000		format(A15,1X, I4,4X, 3(F8.3)) 

!	This gives us the sum of the molecular polarizabilities	Pi(M)
	Ptensor = Ptensor + PolPrime


!	Next, to find the sum of the interaction polarizabilities Pi(I)
	j = j + 1	
		
	endif
	end do
	
	write(12,FMT=14000) 'Ptensor Diagonal', j, Ptensor(1,1), Ptensor(2,2), Ptensor(3,3)
14000		format(A17,1X, I4,4X, 3(F10.3))  
	
end program EulerRotCoor 

!Subroutines	
Subroutine CenterMass (x1,y1,z1, x2,y2,z2, x3,y3,z3, CMass) 
	! input arguments
	implicit none
	!integer, intent(in)	:: j ! each j is a water molecule
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
	real, dimension(3,3) :: Polar(3,3), PolPrime(3,3)
	
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

	PolPrime = InvEuler*polar*Euler
	
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

	
		! write(*,*) 'Molecular Frame X-axis', j 
		! write(*,*) 'X', xaxx, 'Y', xaxy, 'Z', xaxz		
		! write(*,*) 'Molecular Frame Y-axis', j 
		! write(*,*) 'X', yaxx, 'Y', yaxy, 'Z', yaxz		
		! write(*,*) 'Molecular Frame Z-axis', j 
		! write(*,*) 'X', zaxx, 'Y', zaxy, 'Z', zaxz

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
	
!	DummyNorm(:) = [NormX, NormY, NormZ] ! Array constructor
	
	write(*,*) 'Dummy Normal', j 
	write(*,*) 'X', dNormX, 'Y', dNormY, 'Z', dNormZ	

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
	!real ::  MagVec   
	 
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

