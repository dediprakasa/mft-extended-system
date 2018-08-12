    program Percobaan
    
    use global
    
    implicit none
    
    REAL*8,EXTERNAL::fermi,rtbis,integrand
    integer, allocatable::ind(:,:,:)
    INTEGER::i,j,k,ky,kz,l,m,n,lp,mp,np,counter,iw,pw,is
    REAL*8:: time_start,time_finish,time_elapsed
    
        !   Arrays and variables for LAPACK
    complex*16, dimension(2*2)::work
    integer, dimension(2)::ipiv
    integer::info

    call CPU_TIME(time_start)
    
    eps=0.D0
    t1=1.D0
    t2=0.D0
    Nw=501
    wmin=-8.D0
    wmax=8.D0
    nrows=2
    nfilling=nrows/2.D0
    U = 0.d0
    alpha = 0.1d0
    Nk=51



    
    call allocation
    allocate(ind(0:4,0:4,0:4))  

    !   Pemberat integrasi Simpson komposit
    dw=(wmax-wmin)/real(Nw-1,8)
    wfreq(1)=dw/3.D0
    do i=2,Nw-1,2
        wfreq(i)=4.D0*dw/3.D0
    end do
    do i=3,Nw-2,2
        wfreq(i)=2.D0*dw/3.D0
    end do
    wfreq(Nw)=dw/3.D0

    
    !   Membuat Matriks Identitas
    IdMat(:,:)=0.D0
    do i=1,nrows
        IdMat(i,i)=1.D0
    end do
    
    !   Generate Hamiltonian matrix

    H(:,:)=0.D0    
    
    ! Matriks sigma awal
    sigma(:,:,:) = 0.d0

    do j = 1,nrows/2
        n_up(j) = 0.5d0
        n_down(j+nrows/2) = 0.5d0
    end do

    do j = 1,nrows/2
        sigma(j,j,:) = U*n_down(j)
        sigma(j+nrows/2,j+nrows/2,:) = U*n_up(j)
    end do


    ! kkkkkkkkkkkkkkkkkk

    dk=(2.d0*pi/C)/real(Nk-1,8)

    wk(:)=1.d0

    do i=1,Nk
        do j=1,Nk
            do l=1,Nk
                n=(i-1)*(Nk**2)+((j-1)*Nk+l)
                kk(1,n)=(-pi/C)+dk*real(i-1,8)   
                kk(2,n)=(-pi/C)+dk*real(j-1,8)   
                kk(3,n)=(-pi/C)+dk*real(l-1,8)

                if((i==1).or.(i==Nk)) wk(n)=0.5d0*wk(n)
                if((j==1).or.(j==Nk)) wk(n)=0.5d0*wk(n)
                if((l==1).or.(l==Nk)) wk(n)=0.5d0*wk(n)
            enddo
        enddo
    enddo
    normk=sum(wk)


    H(1,2)=0.d0
    H(2,1)=0.d0

    if (U == 0.d0) then
        ism = 1
    else
        ism = 400
    end if

    ! Iterasi
    do is = 1,ism
    
        norm = 0.d0
        print *, 'Iterasi ke ', is
        
        do iw=1,Nw
            sumk(:,:)=(0.d0,0.d0)
            do k=1,Nk
                do ky=1,Nk
                    do kz=1,Nk
                        w(iw)=wmin+dw*(iw-1)
                        H(1,1)=2*t*(cos((kk(1,k)*C))+cos((kk(2,ky)*C))+cos((kk(3,kz)*C)))
                        H(2,2)=H(1,1)       
                        Gdummy(:,:)=(w(iw)+ii*0.05D0)*IdMat(:,:)-H(:,:)-sigma(:,:,iw)
                    end do
                end do
            enddo
                
                !      Matrix Inversion
                call zgetrf(nrows,nrows,Gdummy,nrows,ipiv,info)
                call zgetri(nrows,Gdummy,nrows,ipiv,work,2*nrows,info)

            do k=1,Nk**3
                G(:,:,iw)=Gdummy(:,:)
                sumk(:,:)=sumk(:,:)+G(:,:,iw)*wk(k)
            enddo

                G(:,:,iw)=(1/normk)*sumk(:,:)
                traceup = 0.d0
                tracedown = 0.d0
                
                do j=1,nrows/2
                    traceup = traceup + aimag(G(j,j,iw))
                end do

                do j = (nrows/2)+1,nrows
                    tracedown = tracedown + aimag(G(j,j,iw))
                end do

                do j = 1,nrows
                    PDOS(j,iw)=-(1.D0/pi)*aimag(G(j,j,iw))                      
                end do

                DOSup(iw) = -(1.D0/pi)*(traceup)
                DOSdown(iw) = -(1.D0/pi)*(tracedown)
                DOS(iw)=-(1.D0/pi)*(traceup+tracedown)
                norm=norm+DOS(iw)*wfreq(iw)
            
        end do

        mu=rtbis(integrand,wmin,wmax,tolerance)

        print*, mu


        !   Hitung <n>
        nn(:)=0.D0
        do j=1,nrows  
            do iw = 1,Nw
                nn(j) = nn(j)+PDOS(j,iw)*fermi(mu,T,w(iw))*wfreq(iw)          
            end do  
        end do


        sigmac(:,:)=0.d0
        do iw = 1,Nw
            
            do j = 1,nrows/2
                sigmac(j,j)=U*nn(j+nrows/2)
                sigmac(j+nrows/2,j+nrows/2)=U*nn(j)
            end do

            sigmaf(:,:,iw) = sigmac(:,:)

        end do
        
        
        ! Cek error
        print *, 'error', maxval(abs(sigmaf(:,:,:) - sigma(:,:,:)))
        print *, 'norm',norm
        if (maxval(abs(sigmaf(:,:,:) - sigma(:,:,:))) .lt. 0.008d0) then 
            print *, 'Konvergen'
!           open(unit=111,file='GF.dat',status='unknown')
            open(unit=11,file='DoS_vs_w.dat',status='unknown')
            open(unit=15,file='DoSup_vs_w.dat',status='unknown')
            open(unit=16,file='DoSdown_vs_w.dat',status='unknown')
            open(unit=12,file='mu.dat',status='unknown')
            open(unit=14,file='nnn.dat',status='unknown')
            write(12,*)mu,0.D0
            write(12,*)mu,nrows/2.D0        
                do iw = 1,Nw
                    write(11,*)w(iw),DoS(iw)
                    write(15,*)w(iw),DoSup(iw)
                    write(16,*)w(iw),DOSdown(iw)
 !                   do j = 1,nrows/2
 !                       write(14,*)j,sigma(j+nrows/2,j+nrows/2,iw)/U,sigma(j,j,iw)/U
 !                   end do
                end do
            call Optical_Response

            call CPU_TIME(time_finish)
        
            time_elapsed= time_finish-time_start

            write(*,*) 'All calculations DONE !!!'
            write(*,*)
            write(*,*) 'elapsed time (seconds) =', time_elapsed
            write(*,*) 'elapsed time (minutes) =', time_elapsed/60.0
            write(*,*) 'elapsed time (hours) =', time_elapsed/3600.0
            write(*,*)

            close(11)
            close(15)
            close(16)
            close(12)
            close(14)

            stop
        else
            sigma(:,:,:) = alpha*sigmaf(:,:,:) + (1-alpha)*sigma(:,:,:)
        end if

    
    end do



    deallocate(ind)
    call deallocation



    
    end program
    

    FUNCTION integrand(x)
    
    USE global
    IMPLICIT NONE
    INTEGER :: i
    REAL*8::integrand,x,tot
    REAL*8,EXTERNAL::fermi

    tot=0.D0
    DO i=1,Nw
        tot=tot+DOS(i)*fermi(x,T,w(i))*wfreq(i)
    END DO
    integrand=nfilling-(tot*Real(nrows,8)/norm)
    
    RETURN

    
    END FUNCTION integrand
