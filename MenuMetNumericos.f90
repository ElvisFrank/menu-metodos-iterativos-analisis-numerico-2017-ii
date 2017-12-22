!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			ELVIS FRANK DOMINGUEZ VIDAL						!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!					contenedor de metodos					!
module subrutinas
	!USE, intrinsic :: iso_fortran_env, ONLY: WP => REAL64
	use funciones
	contains
	Subroutine MetBiseccion
	! ---------------------------------------------------
	! METODO DE BISECCION para encontrar una solución
	! de f(x)=0 dada la función continua f en el intervalo
	! [a,b] donde f(a) y f(b) tienen signos opuestos.
	! ---------------------------------------------------

		real::a,b,error,TOL,p,itr,xa,xb,l
		integer i

		i=1;
		!write(*,*)'Ingrese inicio del intervalo, a: '
		!read(*,*)a	!leer a
		a=1
		!write(*,*)'Ingrese fin del intervalo, b: '
		!read(*,*)b	!leer b
		b=2 ![a,b]

		!write(*,*)'Ingrese la tolerancia: '
		!read(*,*)tol	!leer tol
		tol=0.0001
		xa=(b-a)/tol
		xb=2
		IF (fb(a)*fb(b)>0.0) THEN
			PRINT*, 'Es posible que no haya raiz en este intervalo'
		ELSEIF (fb(a)*fb(b)<0.0) THEN
			DO i=1,50
				!print*, 'i: ',i
				p=((a+b)/2.0)
				error=abs((b-a)/b)
				If ((fb(p)==0.0).OR.(error<TOL)) THEN !.or.(o)
					!Mostrando lasa raices encontradas
					Print*, ' La Raiz es = '
					Print'(F15.9)',p
					PRINT*, 'La raiz es',p
					EXIT
				END IF
				!reasignando valores a, b
				IF (fb(a)*fb(p)>0.0) THEN
					a=p
				ELSE IF (fb(a)*fb(p)<0.0) THEN
					b=p
				END IF
			END DO
		END IF

		l=alog10(xa)
		itr=(l)/(alog10(xb)) !iteraciones
		print*,'|----------------------------------------------|'
		print*,'numero de iteraciones'
		print*,'|----------------------------------------------|'
		write(*,2),itr
		write(*,3),itr
		2 format('iteracion',2x,'>',F9.6,//)
		3 format('iteracion',2x,'>',E10.5,//)     
	      
		  
	End Subroutine MetBiseccion
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Subroutine MetNewton
	! ---------------------------------------------------
	! METODO DE NEWTON Para encontrar una solución
	! de f(x) = 0 dada una aproximación inicial p0 :
	! ---------------------------------------------------
		!implicit none <obliga a declarar todas las variables>
		real::tol,po,p,error,er
		integer i,bandera, No
		bandera=0
		No=50
		write(*,*)'elige una tolerancia de error'
		read(*,*)tol	!leer tolerancia
		write(*,*)'elige el punto de comienzo'
		read(*,*)po
		do i=1,No
			print*,i
			p=((po)-(fn(po)/fnd(po)))
			error=((po-p)/(po))	!error
			er=(abs(error)) 	!valor absoluto del error relat.
			if(er<tol) THEN !error abs menor q tolerancia
				bandera=1;
				write(*,*)'el valor de la el valor de la raiz es:',p
				EXIT
			end if
			po=p
		end do
		if(bandera==0)write(*,*)'procedimiento completado sin exito'
	End Subroutine MetNewton
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Subroutine MetNewtonRModf
		real::p0,tol,p1,y,yd,ydd,error
		integer::i,bandera=0, No=50

		!write(*,*)'elige una tolerancia de error'
		!read(*,*)tol	!leer tolerancia
		!write(*,*)'elige el punto de comienzo'
		!read(*,*)p0

		p0=1
		tol=0.001

		do i=1,No
			print*,i
			y=fnm(p0)		!f(p0)
			yd=fnmd(p0)		!f'(p0)
			ydd=fnmdd(p0)	!ff''(p0)
			p1=p0-(y*yd)/((yd**2)-(y*ydd))

			error=abs((p1-p0)/p1)	! |Er|
			if(error<=tol)then
				write(*,*)'la raiz en es: ',p0
				bandera=1;
				EXIT
			end if
			p0=p1
		end do
		if(bandera==0)write(*,*)'procedimiento completado sin exito'
	
	End Subroutine MetNewtonRModf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Subroutine MetSecante
	! ---------------------------------------------------
	! METODO DE LA SECANTE Para encontrar una solución
	! de f(x) = 0 dada una función f en el intervalo
	![a,b] donde f(a) y f(b) tienen signos opuestos:
	! ---------------------------------------------------
		real::a,b,tol,p0,q0,q1,p,error
		integer::i,bandera=0, No=50


		!write(*,*)'elige una tolerancia de error'
		!read(*,*)tol	!leer tolerancia
		!write(*,*)'elige el inicio del intervalo, a:'
		!read(*,*)a
		!write(*,*)'elige el fin del intervalo, b:'
		!read(*,*)b

		a=1.3
		b=1.4
		tol=0.0002
		if(fs(a)>0)then
			p0=b
			q0=fs(a)
			q1=fs(b)
		else if(fs(a)<0) then
			p0=a
			q0=fs(b)
			q1=fs(a)
		end if
		do i=2 ,No			
			if(fs(a)>0)then
				p=p0-(q1/(q1-q0))*(p0-a)
			else if(fs(a)<0) then
				p=p0-(q1/(q0-q1))*(b-p0)
			end if

			error=abs((p-p0)/p)
			if(error<tol)then
				print*,'La solucion aproximada es: ',p
				bandera=1;
				EXIT
			end if
			!redefiniendo p0, q1
			p0=p
			q1=fs(p)
			!write(*,*)'valor de i: ',i
		end do
		if(bandera==0)write(*,*)'procedimiento completado sin exito'

	End Subroutine MetSecante
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Subroutine MetSecanteModf
	! ---------------------------------------------------
	! METODO DE LA SECANTE MODIFICADA Para encontrar
	! una solución de f(x) = 0 dada las aproximaciones
	! iniciales p0 y p1
	! ---------------------------------------------------
		real::p0,p1,tol,q0,q1,p,error
		integer::i,bandera=0, No=50
		p0=1.3
		p1=1.4
		tol=0.0002

		q0=fs(p0)
		q1=fs(p1)

		do i=2 ,No			
			p=p1-(q1/(q1-q0))*(p1-p0)
			error=abs((p-p1)/p)
			if(error<tol)then
				print*,'La solucion aproximada es: ',p
				bandera=1;
				EXIT
			end if
			!redefiniendo p0, q1
			p0=p1
			q0=q1
			p1=p
			q1=fs(p)
			write(*,*)'valor de i: ',i
		end do
		if(bandera==0)write(*,*)'procedimiento completado sin exito'

	End Subroutine MetSecanteModf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Subroutine MetInterpolacionL
		implicit none

	End Subroutine MetInterpolacionL
end module subrutinas

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!					contenedor de funciones					!
module funciones	
	contains

	function fb(x)	!fumcion biseccion
		real,intent(in)::x
		real, parameter::pi=3.1415926536
		real::fb
		!fb=((x*exp(x))-pi) 		!0 2
		fb=((x**3)+(4*(x**2))-10) 	!1 2
		!return
	end function fb

	function fn(x)	!funcion newton
		real,intent(in)::x
		real::fn
		fn=((x**3)+(4*(x**2))-10)
	end function fn

	function fnd(x)	!funcion newton derivada
		real,intent(in)::x
		real::fnd
		fnd=(3*(x**2))+(8*(x))
	end function fnd


	function fnm(x)		!funcion newton 
		real,intent(in)::x
		real::fnm
		fnm=((x**3)+(4*(x**2))-10)
	end function fnm

	function fnmd(x)	!funcion newton derivada
		real,intent(in)::x
		real::fnmd
		fnmd=(3*(x**2))+(8*(x))
	end function fnmd

	function fnmdd(x)	!funcion newton 2da derivada
		real,intent(in)::x
		real::fnmdd
		fnmdd=((6*x)+8)
	end function fnmdd



	function fs(x)		!funcion secante
		real,intent(in)::x
		real::fs
		fs=((x**3)+(4*(x**2))-10)!1.3 1.4
	end function fs

end module funciones

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!						menu principal						!
program mNumerico
	use subrutinas
	use funciones

	implicit none
	integer:: opc
	character:: cont

	156 continue
	print*, '|----------------------------------------------|'
	print*, '			METODOS ITERATIVOS'
	print*, '|----------------------------------------------|'
	print*, '	1 -METODO BISECCION'
	print*, '	2 -METODO NEWTON'
	print*, '	3 -METODO NEWTON RAPSON MODIFICADO'
	print*, '	4 -METODO SECANTE'
	print*, '	5 -METODO SECANTE MODIFICADO'
	print*, '	6 -METODO INTERPOLACION LINEAL'
	print*, '|----------------------------------------------|'
	read*, opc

	select case(opc)
		case(1)
			call MetBiseccion
		case(2)
			call MetNewton
		case(3)
			call MetNewtonRModf
		case(4)
			call MetSecante
		case(5)
			call MetSecanteModf
		case(6)
			call MetInterpolacionL
		case default
		print*, '	¡¡¡La opcion que ingreso NO ES VALIDA!!!'

		goto 156
	End select
	print*, '|----------------------------------------------|'
	print*, '	DESEA OTRA VEZ EJECUTAR AL PROGRAMA S/N'
	Print*, ''
	read*, cont
	if (cont=='s')then
		goto 156
	End if

end program mNumerico