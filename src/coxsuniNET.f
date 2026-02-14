      subroutine coxsuninet(parm,nobs,nvars,x,y,d,g,w,jd,vp,ne,n
     *x,nlam,flmin,ulam,thr, maxit,isd,lmu,ca,ia,nin,dev0,dev,a
     *lm,nlp,jerr,iglampos, olampos, alpha)
      implicit double precision(a-h,o-z)
      double precision x(nobs,nvars),y(nobs),d(nobs),g(nobs),w(
     *nobs),vp(nvars),ulam(nlam)
      double precision ca(nx,nlam),dev(nlam),alm(nlam)
      integer jd(*),ia(nx),nin(nlam)
      double precision, dimension (:), allocatable :: xs,ww,vq
      integer, dimension (:), allocatable :: ju
      if(maxval(vp) .gt. 0.0)goto 10021
      jerr=10000
      return
10021 continue
      allocate(ww(1:nobs),stat=jerr)
      if(jerr.ne.0) return
      allocate(ju(1:nvars),stat=jerr)
      if(jerr.ne.0) return
      allocate(vq(1:nvars),stat=jerr)
      if(jerr.ne.0) return
      if(isd .le. 0)goto 10041
      allocate(xs(1:nvars),stat=jerr)
      if(jerr.ne.0) return
10041 continue
      call chkvars(nobs,nvars,x,ju)
      if(jd(1).gt.0) ju(jd(2:(jd(1)+1)))=0
      if(maxval(ju) .gt. 0)goto 10061
      jerr=7777
      return
10061 continue
      vq=max(0d0,vp)
      vq=vq*nvars/sum(vq)
      ww=max(0d0,w)
      sw=sum(ww)
      if(sw .gt. 0.0)goto 10081
      jerr=9999
      return
10081 continue
      ww=ww/sw
      call cstandard(nobs,nvars,x,ww,ju,isd,xs)
      if(isd .le. 0)goto 10101
      do 10111 j=1,nvars
10111 continue
      continue
10101 continue
      call coxsuninet1(parm,nobs,nvars,x,y,d,g,ww,ju,vq,ne,nx,n
     *lam,flmin,ulam, thr,  isd,maxit,lmu,ca,ia,nin,dev0,dev,al
     *m,nlp,jerr,iglampos, olampos, alpha)
      if(jerr.gt.0) return
      dev0=2.0*sw*dev0
      if(isd .le. 0)goto 10131
      do 10141 k=1,lmu
      nk=nin(k)
      ca(1:nk,k)=ca(1:nk,k)/xs(ia(1:nk))
10141 continue
      continue
10131 continue
      deallocate(ww,ju,vq)
      if(isd.gt.0) deallocate(xs)
      return
      end

      subroutine coxsuninet1(parm,nobs,nvars,x,y,d,g,q,ju,vp,n
     *e,nx,nlam,flmin,ulam,cthri,  isd,maxit,lmu,ao,m,kin,dev0
     *,dev,alm,nlp,jerr,iglampos, olampos, alpha)
      implicit double precision(a-h,o-z)
      double precision x(nobs,nvars),y(nobs),q(nobs),d(nobs)
      double precision g(nobs),vp(nvars),ulam(nlam)
      double precision ao(nx,nlam),dev(nlam),alm(nlam)
      integer ju(nvars),m(nx),kin(nlam)
      double precision, dimension (:), allocatable :: w,dk,v,x
     *s,wr
      double precision, dimension (:), allocatable :: a,as,f,dq
      double precision, dimension (:), allocatable :: e,uu,ga
      integer, dimension (:), allocatable :: jp,kp,mm,ixx
      call get_int_parms(sml,eps,big,mnlam,devmax,pmin,exmx,
     *itrace)
      isd = isd*1
      sml=sml*100.0
      devmax=devmax*0.99/0.999
      allocate(e(1:nobs),stat=jerr)
      if(jerr.ne.0) return
      allocate(uu(1:nobs),stat=jerr)
      if(jerr.ne.0) return
      allocate(f(1:nobs),stat=jerr)
      if(jerr.ne.0) return
      allocate(w(1:nobs),stat=jerr)
      if(jerr.ne.0) return
      allocate(v(1:nvars),stat=jerr)
      if(jerr.ne.0) return
      allocate(a(1:nvars),stat=jerr)
      if(jerr.ne.0) return
      allocate(as(1:nvars),stat=jerr)
      if(jerr.ne.0) return
      allocate(xs(1:nvars),stat=jerr)
      if(jerr.ne.0) return
      allocate(ga(1:nvars),stat=jerr)
      if(jerr.ne.0) return
      allocate(ixx(1:nvars),stat=jerr)
      if(jerr.ne.0) return
      allocate(jp(1:nobs),stat=jerr)
      if(jerr.ne.0) return
      allocate(kp(1:nobs),stat=jerr)
      if(jerr.ne.0) return
      allocate(dk(1:nobs),stat=jerr)
      if(jerr.ne.0) return
      allocate(wr(1:nobs),stat=jerr)
      if(jerr.ne.0) return
      allocate(dq(1:nobs),stat=jerr)
      if(jerr.ne.0) return
      allocate(mm(1:nvars),stat=jerr)
      if(jerr.ne.0) return
      call groups(nobs,y,d,q,nk,kp,jp,t0,jerr)
      if(jerr.ne.0) go to 10180
      oma=1.0-parm
      nlm=0
      ixx=0
      altemp=0.0
      dq=d*q
      call died(nobs,nk,dq,kp,jp,dk)
      a=0.0
      f(1)=0.0
      fmax=log(huge(f(1))*0.1)
      if(nonzero(nobs,g) .eq. 0)goto 10201
      f=g-dot_product(q,g)
      e=q*exp(sign(min(abs(f),fmax),f))
      goto 10211
10201 continue
      f=0.0
      e=q
10211 continue
      continue
      r0=risk(nobs,nvars,nk,dq,dk,f,e,kp,jp,uu)
      rr=-(dot_product(dk(1:nk),log(dk(1:nk)))+r0)
      dev0=rr
      do 10221 i=1,nobs
      if((y(i) .ge. t0) .and. (q(i) .gt. 0.0))goto 10241
      w(i)=0.0
      wr(i)=w(i)
10241 continue
10221 continue
      continue
      call outer(nobs,nk,dq,dk,kp,jp,e,wr,w,jerr,uu)
      if(jerr.ne.0) go to 10180
      alf=1.0
      if(flmin .ge. 1.0)goto 10261
      eqs=max(eps,flmin)
      alf=eqs**(1.0/(nlam-1))
10261 continue
      m=0
      mm=0
      nlp=0
      nin=nlp
      mnl=min(mnlam,nlam)
      as=0.0
      cthr=cthri*dev0
      do 10271 j=1,nvars
      if(ju(j).eq.0)goto 10271
      ga(j)=abs(dot_product(wr,x(:,j)))
10271 continue
      continue
      do 10281 ilm=1,nlam
      ! if(itrace.ne.0) call setpb(ilm-1)
      al0=altemp
      if(flmin .lt. 1.0)goto 10301
      altemp=ulam(ilm)
      al = altemp * (1.0 - alpha)
      goto 10291
10301 if(ilm .le. 2)goto 10311
      altemp=altemp*alf
      al = altemp * (1.0 - alpha)
      goto 10291
10311 if(ilm .ne. 1)goto 10321
      altemp=big
      al = big
      goto 10331
10321 continue
      al0=0.0
      do 10341 j=1,nvars
      if(ju(j).eq.0)goto 10341
      if(vp(j).gt.0.0) al0=max(al0,ga(j)/vp(j))
10341 continue
      continue
      al0=al0/max(parm,1.0d-3)
      altemp=alf*al0
      al = altemp * (1.0 - alpha)
10331 continue
10291 continue
      if (iglampos .eq. 1) olampos = altemp * alpha * 2.0
      if (iglampos .eq. 1) al = al * 2.0
      omal=oma * altemp
      tlam=parm * (2.0*altemp-al0)
      do 10351 k=1,nvars
      if(ixx(k).eq.1)goto 10351
      if(ju(k).eq.0)goto 10351
      if(ga(k).gt.tlam*vp(k)) ixx(k)=1
10351 continue
      continue
10360 continue
      continue
10371 continue
      if(nlp .le. maxit)goto 10391
      jerr=-ilm
      return
10391 continue
      if(nin.gt.0) as(m(1:nin))=a(m(1:nin))
      call vars(nobs,nvars,x,w,ixx,v)
      continue
10401 continue
      nlp=nlp+1
      dli=0.0
      do 10411 j=1,nvars
      if(ixx(j).eq.0)goto 10411
      u=a(j)*v(j)+dot_product(wr,x(:,j))
      at = 0.0
      if(u .gt. vp(j)*al) at = (u - vp(j)*al)/(v(j)+vp(j)*omal)
      if(u .lt. - olampos) at = (u + olampos)/(v(j)+vp(j)*omal)
!if(abs(u) .gt. vp(j)*al)goto 10431
!at=0.0
!goto 10441
! 10431 continue
!at=sign(abs(u)-vp(j)*al,u)/(v(j)+vp(j)*omal)
! 10441 continue     
      continue
      if(at .eq. a(j))goto 10461
      del=at-a(j)
      a(j)=at
      dli=max(dli,v(j)*del**2)
      wr=wr-del*w*x(:,j)
      f=f+del*x(:,j)
      if(mm(j) .ne. 0)goto 10481
      nin=nin+1
      if(nin.gt.nx)goto 10412
      mm(j)=nin
      m(nin)=j
10481 continue
10461 continue
10411 continue
10412 continue
      if(nin.gt.nx)goto 10402
      if(dli.lt.cthr)goto 10402
      if(nlp .le. maxit)goto 10501
      jerr=-ilm
      return
10501 continue
      continue
10511 continue
      nlp=nlp+1
      dli=0.0
      do 10521 l=1,nin
      j=m(l)
      u=a(j)*v(j)+dot_product(wr,x(:,j))
      at = 0.0
      if(u .gt. vp(j)*al) at = (u - vp(j)*al)/(v(j)+vp(j)*omal)
      if(u .lt. - olampos) at = (u + olampos)/(v(j)+vp(j)*omal)
!if(abs(u) .gt. vp(j)*al)goto 10541      
!at=0.0
!goto 10551   
! 10541 continue     
!at=max(cl(1,j),min(cl(2,j),sign(abs(u)-vp(j)*al,u)/  (v(j)+vp(j)*o
!      *mal)))
! 10551 continue     
      continue
      if(at .eq. a(j))goto 10571
      del=at-a(j)
      a(j)=at
      dli=max(dli,v(j)*del**2)
      wr=wr-del*w*x(:,j)
      f=f+del*x(:,j)
10571 continue
10521 continue
      continue
      if(dli.lt.cthr)goto 10512
      if(nlp .le. maxit)goto 10591
      jerr=-ilm
      return
10591 continue
      goto 10511
10512 continue
      goto 10401
10402 continue
      if(nin.gt.nx)goto 10372
      e=q*exp(sign(min(abs(f),fmax),f))
      call outer(nobs,nk,dq,dk,kp,jp,e,wr,w,jerr,uu)
      if(jerr .eq. 0)goto 10611
      jerr=jerr-ilm
      go to 10180
10611 continue
      ix=0
      do 10621 j=1,nin
      k=m(j)
      if(v(k)*(a(k)-as(k))**2.lt.cthr)goto 10621
      ix=1
      goto 10622
10621 continue
10622 continue
      if(ix .ne. 0)goto 10641
      do 10651 k=1,nvars
      if(ixx(k).eq.1)goto 10651
      if(ju(k).eq.0)goto 10651
      ga(k)=abs(dot_product(wr,x(:,k)))
      if(ga(k) .le. al*vp(k))goto 10671
      ixx(k)=1
      ix=1
10671 continue
10651 continue
      continue
      if(ix.eq.1) go to 10360
      goto 10372
10641 continue
      goto 10371
10372 continue
      if(nin .le. nx)goto 10691
      jerr=-10000-ilm
      goto 10282
10691 continue
      if(nin.gt.0) ao(1:nin,ilm)=a(m(1:nin))
      kin(ilm)=nin
      alm(ilm)=altemp
      lmu=ilm
      dev(ilm)=(risk(nobs,nvars,nk,dq,dk,f,e,kp,jp,uu)-r0)/rr
      if(ilm.lt.mnl)goto 10281
      if(flmin.ge.1.0)goto 10281
      me=0
      do 10701 j=1,nin
      if(ao(j,ilm).ne.0.0) me=me+1
10701 continue
      continue
      if(me.gt.ne)goto 10282
      if((dev(ilm)-dev(ilm-mnl+1))/dev(ilm).lt.sml)goto 10282
      if(dev(ilm).gt.devmax)goto 10282
10281 continue
10282 continue
      g=f
10180 continue
      deallocate(e,uu,w,dk,v,xs,f,wr,a,as,jp,kp,dq,mm,ga,ixx)
      return
      end