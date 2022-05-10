program CMassTest  

	implicit none 

	real H1Wx, H1Wy, H1Wz 
	real H2Wx, H2Wy, H2Wz
	real OWx, OWy, OWz 

	real dNormX, dNormY, dNormZ
	real O1HvecX, O1HvecY, O1HvecZ !dipole vector/bond lengths along the H1W to O bond
	real O2HvecX, O2HvecY, O2HvecZ !dipole vector/bond lengths along the H2W to O bond

	real xdirx, xdiry, xdirz, xaxx, xaxy, xaxz !Molecular Frame x-axis
	real yaxx, yaxy, yaxz !Molecular Frame y-axis
	real zdirx, zdiry, zdirz, zaxx, zaxy, zaxz !Molecular Frame z-axis
	real magHH, normXproj, HHtheta, normO1H, normO2H
	real HHx, HHy, HHz
	real xout,yout,zout
	real Zlabx, Zlaby, Zlabz, Ylabx, Ylaby, Ylabz, Xlabx, Xlaby, Xlabz
	real phi, theta, psi
	real, dimension(3,3) :: matinv3
	real alpha(3), CMass(3)
	real, dimension(3,3) :: Euler, InvEuler, PolPrime
	integer j

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

!This is the basics of the molecular frame subroutine to make sure it works
!molframe(i, j, H1Wx, H1Wy, H1Wz, H2Wx, H2Wy, H2Wz, NormFact, 
!& dNormX, dNormY, dNormZ, O1HvecX, O1HvecY, O1HvecZ, O2HvecX, O2HvecY, O2HvecZ)
	! atom postions from the water
	H1Wx =26.663
	H1Wy = 14.654
	H1Wz = 39.047 
	H2Wx = 26.225
	H2Wy = 16.226
	H2Wz = 39.124
	OWx = 26.994
	OWy = 15.598
	OWz = 39.010
	! Vectors from the hydrogens to the oxygen
	O1HvecX = -0.331 
	O1HvecY = -0.944 
	O1HvecZ = 0.037
	O2HvecX = -0.769
	O2HvecY = 0.628 
	O2HvecZ = 0.114
	! cross product of the OH vectors for the y-axis
	dNormX = -0.131 
	dNormY = 0.009
	dNormZ = -0.933
	
	! iteration counter that is also the molecule number
	j = 1

	! The oxygen atom is set as the origin for the molecular frame.
	
	! x-axis
	! This is calculated through finding the line that bisects the triangle
	! formed by the H-O-H bond through the oxygen in the molecular frame. 
	! need O1HVec and O2Hvec and the VecMag
	
	call xaxis(O1HvecX, O1HvecY, O1HvecZ, O2HvecX, O2HvecY, O2HvecZ, xout,yout,zout)
		xaxx = xout
		xaxy = yout
		xaxz = zout

	! z-axis 
	! This is calculated by finding the vector between the hydrogens and
	! translating the z variable to the oxygen position in the molecular frame. 
	call zaxis(H1Wx, H1Wy, H1Wz, H2Wx, H2Wy, H2Wz, xaxz, xout, yout, zout)
		zaxx = xout
		zaxy = yout
		zaxz = zout

	! y-axis
	! The y-axis is orthogonal to the z-axis and x-axis such that the Oxygen
	! lone pairs in the classic vespr model of water are in the XY plain 
	! with the C2v symmetry preserved.
		yaxx = dNormX
		yaxy = dNormY 
		yaxz = dNormZ
	
	
			
	write(*,*) 'Molecular Frame X-axis', j 
	write(*,*) 'X', xaxx, 'Y', xaxy, 'Z', xaxz		
	write(*,*) 'Molecular Frame Y-axis', j 
	write(*,*) 'X', yaxx, 'Y', yaxy, 'Z', yaxz		
	write(*,*) 'Molecular Frame Z-axis', j 
	write(*,*) 'X', zaxx, 'Y', zaxy, 'Z', zaxz
	
	! My own coding of the Euler angles 
	! We save this part in the output file
	call EulerAng(xaxx, xaxy, xaxz, Zlabx, Zlaby, Zlabz, zaxx, zaxy, zaxz, &
					Xlabx, Xlaby, Xlabz, phi, theta, psi) 
	write (*,*) 'phi ', phi, 'theta', theta, 'psi', psi
	
	! My own coding of the Euler Rotation matrix and its inverse.
	! This bit gets discarded for output file as unnecessary to keep.
	call EulerRot(phi, theta, psi, Euler, InvEuler)
	write (*,*) 'Euler Rotation Matrix', Euler
	write (*,*) 'Inverse of Euler Rotation matrix', InvEuler
	
	! My own coding of PolarRot = InvEuler*polar*Euler. 
	! This is the bit we save in the output file because it sums later.
	call PolarRot (Euler, InvEuler, alpha, PolPrime)
	write (*,*) 'InvEuler*polar*Euler = PolPrime', PolPrime
	
	call CenterMass(H1Wx, H1Wy, H1Wz, H2Wx, H2Wy, H2Wz, OWx, OWy, OWz , CMass)
	write (*,*) 'Center of Mass', CMass
	
end program CMassTest 

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
			zdirz = HHz + xaxz
			magHH = VecMag(zdirx, zdiry, zdirz)
	
			xout = HHx / magHH
			yout = HHy / magHH
			zout = HHz / magHH
			
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


!Functions
 ! pure function matinv3(A) result(InvA)
    ! implicit none
	! !! Performs a direct calculation of the inverse of a 3×3 matrix.
    ! real, intent(in) :: A(3,3)   !! Matrix
    ! real             :: InvA(3,3)   !! Inverse matrix
    ! real             :: detinv

    ! ! Calculate the inverse determinant of the matrix
    ! detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
              ! - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
              ! + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

    ! ! Calculate the inverse of the matrix
    ! InvA(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    ! InvA(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
    ! InvA(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    ! InvA(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
    ! InvA(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    ! InvA(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
    ! InvA(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    ! InvA(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
    ! InvA(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
 ! end function
