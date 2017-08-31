 program main
	 !this program is a 2D linear case, it is useful for a reaction step solution as it allows smooth of the function.
 IMPLICIT none
 integer,parameter::jx=200
 integer,parameter::jy=200
 double precision u(-2:jx+2,-2:jy+2)
 double precision ax(-2:jx+2,-2:jy+2)
 double precision ay(-2:jx+2,-2:jy+2) 
 double precision x(-2:jy+2)
 double precision y(-2:jy+2)
 
 double precision px(-2:jx+2,-2:jy+2)        !位置标记函数
 double precision py(-2:jx+2,-2:jy+2)
 double precision px_tmp(-2:jx+2,-2:jy+2)
 double precision py_tmp(-2:jx+2,-2:jy+2)
 double precision V(-2:jx+2,-2:jy+2)
 double precision V_tmp(-2:jx+2,-2:jy+2)
 double precision cx(-2:jx+2,-2:jy+2) 
 double precision cy(-2:jx+2,-2:jy+2) 
 
 double precision Px_init(0:2*jx,0:2*jy)
 double precision Py_init(0:2*jx,0:2*jy)
 integer NnX_init(0:2*jx,0:2*jy)
 integer NNy_init(0:2*jx,0:2*jy)
 double precision V_init(0:2*jx,-jy:2*jy) 
 double precision cx_init(0:2*jx,0:2*jy)
 double precision cy_init(0:2*jx,0:2*jy)
 integer mark(0:jx,0:jy)

 double precision ax_tmp,ay_tmp
 double precision D_W,U_W,L_W,R_W,qq,xx,yy,qq1
 double precision dx,dt,tt,t,CFL,pi,dy,ratio,ss,d,radius
 double precision dist,x0,y0,amax,dlx,dly
 integer i,j,n1,aa,it,ij,ii,nx_init,ny_init,n_new,ix,iy
integer drx,dry,dpx,dpy,dsx,dsy,kind_prblm
double precision L,R

!    n=100
   ss=1.d-6 
    pi=6*dasin(0.5d0)
!    tt=0.6d0
!
kind_prblm=2

!**************************************************************!!************初始条件***************
select case(kind_prblm)
case(1)
	dlx=2.d0
	dly=2.d0
	dx=dlx/jx
	dy=dly/jy

	CFL=0.2d0
	tt=1.11
	!tt=0.5
	u=100
  amax=0
   do i=-2,jx+2
       do j=-2,jy+2
            x(i)=dx*i
            y(j)=dy*j
 
			ax(i,j)=1
			ay(i,j)=1
 amax=1
			if(x(i).ge.0.1.and.x(i).le.0.5.and.y(j).le.0.5.and.y(j).ge.0.1)then
				u(i,j)=1
			else 
				u(i,j)=0
			endif
		enddo
enddo
!rotation Circle
!
!    ss=1.d-10
!    dx=1.d0/jx
!    dy=1.d0/jy
!
!    tt=1.
!
!    CFL=0.02
!    dt=CFL*min(dx,dy)/pi
!    do i=-2,jx+2
!        x0=i*dx
!		do j=-2,jy+2
!            y0=j*dy
!			r=sqrt((x0-0.5)**2+(y0-0.5)**2)
!			if(r.ge.0.4)then
!               u(i,j)=0
!		   else if(y0.gt.0.3999.and.y0.lt.0.5999.and.x0.ge.0.5)then
!               u(i,j)=0
!			else
!			u(i,j)=1
!		   endif
!    ax(i,j)=-pi*(y0-0.5)
!	ay(i,j)=pi*(x0-0.5)
!        enddo
!    enddo
!
!


!************************************************************


!    
!    do i=-2,jx+2
!        x=i*dx
!        do j=-2,jy+2
!            y=j*dy
!            if(x.ge.0.2.and.x.le.0.5.and.y.ge.0.2.and.y.le.0.5)then
!               u(i,j)=1
!            else
!               u(i,j)=0
!            endif
!            
!        enddo
!    enddo
!   !
!  

   !ax=1
   !ay=1
   !
   !jx=250
   !jy=150
  ! x0=0.5d0
  ! y0=0.3d0
  ! 
  ! dlx=1.0d0
  ! dly=1.0d0
  ! dx=dlx/jx
  ! dy=dly/jy
  ! radius=0.2d0
  ! CFL=0.6d0
  ! tt=1.0d0
 !***********************************  stretching Circle
 case(2)
   x0=pi/2
   y0=0.7d0
   !
   !dlx=1.0d0
   !dly=1.0d0
   !dx=dlx/jx
   !dy=dly/jy
   !radius=0.2d0
   !CFL=0.1d0
   !tt=1.0d0
   
   dlx=pi
   dly=pi
   dx=dlx/jx
   dy=dly/jy
   radius=pi/5
   CFL=0.06d0
   tt=pi
   
   amax=0
   do i=-2,jx+2
       do j=-2,jy+2
            x(i)=dx*i
            y(j)=dy*j
            
            ax(i,j)=dcos((x(i)))*dsin((y(j))) 
            ay(i,j)=-dsin((x(i)))*dcos((y(j))) 
            d=dist(x(i),y(j),x0,y0)
        !if (abs(ay(i,j).ge.amax) then
        if(abs(ay(i,j)).ge.amax) then
            amax=abs(ay(i,j))
       endif
       
            if(d.lt.radius) then 
                u(i,j)=1.0d0 
            else 
                u(i,j)=0.0d0 
            endif
       enddo
   enddo
   end select 
   dt=CFL*min(dx,dy)/amax
     
 !                               ************init***************
   nx_init=jx
   ny_init=jy
   n_new=0
   
   do i=0,jx
       do j=0,jy
        v_init(i,j)=u(i,j)
        px_init(i,j)=0
        py_init(i,j)=0
        nnx_init(i,j)=i
        nny_init(i,j)=j
        px(i,j)=0
        py(i,j)=0
		v(i,j)=u(i,j)
       enddo
   enddo
   
 !  do ii=1,2 
!!***************************!*******
    t=0
    ii=1
    do it=1,1000000                             

!******************time_solve ****begin
        if(t>=tt)then
            exit
        else if(t+dt>tt)then
            dt=tt-t
        endif   
            t=t+dt
        write(*,*) t
       
       mark=0 
 !***********************first_step**********************
 !*************************c_x&c_y***********************
!	if(ii==1)then
    do i=0,jx
        do j=0,jy
            xx=nnx_init(i,j)*dx+px_init(i,j)*dx
            yy=nny_init(i,j)*dy+py_init(i,j)*dy
            cx_init(i,j)=dcos(xx-pi/2)*dsin(yy-pi/2) 
            cy_init(i,j)=-dsin(xx-pi/2)*dcos(yy-pi/2)
		  ! cx_init(i,j)=-pi*(yy-0.5)
		  ! cy_init(i,j)=pi*(xx-0.5)
!		  cx_init(i,j)=1
!		  cy_init(i,j)=1
        enddo
    enddo   
        do i=-1,jx+1
            do j=-1,jy+1
            xx=i*dx+px(i,j)*dx
            yy=j*dy+py(i,j)*dy
            cx(i,j)=dcos(xx-pi/2)*dsin((yy-pi/2)) 
            cy(i,j)=-dsin((xx-pi/2))*dcos(yy-pi/2)
!j			cx(i,j)=1
!			cy(i,j)=1
		  ! cx(i,j)=-pi*(yy-0.5)
		  ! cy(i,j)=pi*(xx-0.5)
			enddo
        enddo
       
       do i=-1,jx+1
            j=-1
         !   if (mark(i,j)==0)then
                mark(i,j)=1
                V(i,j)=0
                px(i,-1)=py(jx+1,jy-i)
                py(i,-1)=px(jx+1,jy-i)
         !   endif
			enddo
        
        do i=-1,jx+1
            j=jy+1
         !j   if (mark(i,j)==0)then
                mark(i,j)=1
                V(i,j)=0
                px(i,jy+1)=py(1,jy-i)
                py(i,jy+1)=px(1,jy-i)
		!	endif
        enddo
        do j=-1,jy+1
            i=-1
            if (mark(i,j)==0)then
                mark(i,j)=1
                V(i,j)=0
                px(-1,j)=py(jx-j,jy-1)
                py(-1,j)=px(jx-j,jy-1)
			endif
        enddo
        do j=-1,jy+1
            i=jx+1
        !    if (mark(i,j)==0)then
                mark(i,j)=1
                V(i,j)=0
                px(i,j)=py(jx-j,1)
                py(i,j)=px(jx-j,1)
        !    endif
        enddo
 
 !******************second_step**************************
 !*******************refresh nn&p************************
        do i=0,jx
            do j=0,jy
                qq=px_init(i,j)+cx_init(i,j)*dt/dx
                nnx_init(i,j)=nnx_init(i,j)+floor(qq+0.5d0)
                px_init(i,j)=qq-floor(qq+0.5d0)
                qq=py_init(i,j)+cy_init(i,j)*dt/dy
                nny_init(i,j)=nny_init(i,j)+floor(qq+0.5d0)
                py_init(i,j)=qq-floor(qq+0.5d0)
            enddo
        enddo
        
        do i=-1,jx+1
            do j=-1,jy+1
            qq=px(i,j)+cx(i,j)*dt/dx
            ix=i+floor(qq+0.5d0)
            qq1=py(i,j)+cy(i,j)*dt/dy
            iy=j+floor(qq1+0.5d0)
                if (ix.ge.0.and.ix.le.jx.and.iy.ge.0.and.iy.le.jy)then
				!	if(mark(ix,iy)==0)then
                    px_tmp(ix,iy)=qq-floor(qq+0.5d0)
                    py_tmp(ix,iy)=qq1-floor(qq1+0.5d0)
                    mark(ix,iy)=1
                    V_tmp(ix,iy)=V(i,j)
			endif
            enddo
        enddo
        
       do i=0,jx
            do j=0,jy
				ix=nnx_init(i,j)
				iy=NNy_init(i,j)
                if(ix.ge.0.and.ix.le.jx.and.iy.ge.0.and.iy.le.jy)then
				!	if(mark(ix,iy).le.1)then
				
                    mark(ix,iy)=2
                    V_tmp(ix,iy)=V_init(i,j)
                    px_tmp(ix,iy)=px_init(i,j)
                    py_tmp(ix,iy)=py_init(i,j)
           endif!
            enddo
       enddo 
      !************************************************************boundary************************************   
     !********************************************************************************   
    ! Solve the interface line 
	!Step 1 find the interface cell and mark them
!	do i=0,jx
!		do j=0,jy
!			if (v(i,j)==1)then
!					if(V(i+ii,j)==0.or.V(i,j+ii)==0)then
!						mark_I(i,j)=1
!						continue
!				enddo



	 !********************************solve blank point******************************
 
    do i=0,jx
        do j=0,jy
            if(mark(i,j)==0)then
				 qq=-cx(i,j)*dt/dx 
				 ix=i+floor(qq+0.5d0)
				 qq=-cy(i,j)*dt/dy 
				 iy=j+floor(qq+0.5d0)
				 v_tmp(i,j)=v(ix,iy)
				 px_tmp(i,j)=0
				 py_tmp(i,j)=0
				endif
    enddo
 enddo
 px=px_tmp
 py=py_tmp
 v=v_tmp
enddo       

    open(1,file='result.plt',status='unknown')
    write(1,*)'TITLE     = "Dataset"'
    write(1,*)'VARIABLES = "x" "y""u"  ZONE T="Zone 1"'
   ! write(1,*)'I=',jx+1,'J=',jy+1,'K=',1,'ZONETYPE=Ordered'
    write(1,*)'DATAPACKING=POINT'

    do 80 i=0,jx
    do 80 j=0,jy

    write(1,*)(i+px(i,j))*dx,(j+py(i,j))*dy,v(i,j)
    80    continue
    close(1)
   read(*,*)it
    end program
    
    
       double precision function dist(x1,y1,x2,y2)  
      implicit double precision(a-h,o-z) 
      dist=dsqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)) 
      end 
