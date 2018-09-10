	MODULE global
	
!	Mathematical and physical constants
	COMPLEX*16,PARAMETER::ii=CMPLX(0.D0,1.D0)
	REAL*8,PARAMETER::PI=2.D0*ASIN(1.D0)
	REAL*8,PARAMETER::tolerance=1.D-5, T=0.D0, C=0.5d-10

	
!	Variables
	REAL*8::eps,t1,t2,t3,wmin,wmax,dw,dws,nfilling,norm,mu,delta_sigma,U,alpha,traceup,tracedown,start,finish,dk,normk
	INTEGER::Nk,nrows,Nw
		
	
!	Matrices
	COMPLEX*16,ALLOCATABLE::G(:,:,:), Gdummy(:,:), Ginv(:,:,:), Gmf_inv(:,:,:), Gloc_inv(:,:,:), sumk(:,:)
	REAL*8,ALLOCATABLE::IdMat(:,:),H(:,:),w(:),DoS(:),PDOS(:,:), sigmad(:,:),wfreq(:),sigma(:,:,:),sigmac(:,:)
	REAL*8,ALLOCATABLE::sigma_down(:,:,:),sigmaf(:,:,:),sigma_up(:,:,:),cek_sigma(:,:,:),nn(:), n_up(:), n_down(:)
	REAL*8,ALLOCATABLE::DOSup(:),DOSdown(:), kk(:,:), wk(:)

		
	END MODULE global

