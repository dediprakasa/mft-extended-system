    subroutine allocation
    
    use global
    implicit none

    allocate(Gdummy(nrows,nrows))
    allocate(G(nrows,nrows,Nw))
    allocate(IdMat(nrows,nrows))
    allocate(H(nrows,nrows))
    allocate(w(Nw))
    allocate(DoS(Nw))
    allocate(PDOS(nrows,Nw))
    allocate(sigmad(nrows,nrows))
    allocate(wfreq(Nw))
    allocate(sigma(nrows,nrows,Nw))
    allocate(sigmac(nrows,nrows))
    allocate(sigmaf(nrows,nrows,Nw))
    allocate(sigma_up((nrows/2)+1,(nrows/2)+1,Nw))
    allocate(sigma_down((nrows/2)+1,(nrows/2)+1,Nw))
    allocate(Ginv(nrows,nrows,Nw))
    allocate(Gmf_inv(nrows,nrows,Nw))
    allocate(Gloc_inv(nrows,nrows,Nw))
    allocate(cek_sigma(nrows,nrows,Nw))
    allocate(nn(nrows))
    allocate(n_up((nrows/2)+1))
    allocate(n_down((nrows/2)+1))
    allocate(DOSup(Nw))
    allocate(DOSdown(Nw))
    allocate(kk(3,Nk**3))
    allocate(wk(Nk**3))
    allocate(sumk(nrows,nrows))

    return
    
    end subroutine allocation


    subroutine deallocation
    
    use global
    implicit none

    deallocate(Gdummy)
    deallocate(G)
    deallocate(IdMat)
    deallocate(H)
    deallocate(w)
    deallocate(DoS)
    deallocate(PDOS)
    deallocate(sigmad)
    deallocate(wfreq)
    deallocate(sigma)
    deallocate(sigmac)
    deallocate(sigmaf)
    deallocate(sigma_up)
    deallocate(sigma_down)
    deallocate(Ginv)
    deallocate(Gloc_inv)
    deallocate(Gmf_inv)
    deallocate(cek_sigma)
    deallocate(nn)
    deallocate(n_up)
    deallocate(n_down)
    deallocate(DOSup)
    deallocate(DOSdown)
    deallocate(wk)
    deallocate(kk)
    deallocate(sumk)


    return
    
    end subroutine deallocation
