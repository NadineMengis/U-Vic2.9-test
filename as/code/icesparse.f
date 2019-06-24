! source file: /net/mare/home1/eby/as/ism/icesparse.F
c===================
c=====

c-----------------------------------------------------------------------

c=================
c=====

c-----------------------------------------------------------------------

c====================
c=====

c-----------------------------------------------------------------------

c===================
c=====

c-----------------------------------------------------------------------

c================
c=====

c-----------------------------------------------------------------------

      subroutine gaussdo (nuvtot, rhs, vec, info)

c     Sets up to do Gaussian elimination. Always called (regardless
c     of defines) if nuvtot <= nuvsmall, to avoid sparse-solution
c     problems for nuvtot <= ~2.

      use comicegrid
      use comicesparse

      dimension rhs(nuvmax), vec(nuvmax)

c     If not defined GAUSS, gaussdo only called for nuvtot <= nuvsmall

      dimension a(nuvsmall,nuvsmall), asav(nuvsmall,nuvsmall)
      nuvgauss=nuvsmall

      info = 0

      call zero (a, nuvgauss*nuvgauss)

      do i=1,nuvtot
        a(i,i) = elspa(i)
      enddo
      do i=1,nuvtot
        ka = ijspa(i)
        kb = ijspa(i+1) -1
        if (kb.ge.ka) then
          do k=ka,kb
            a(i,ijspa(k)) = elspa(k)
          enddo
        endif
      enddo

c       print out in 2-D matrix form

c       solve for vec by Gaussian elimination

      do i=1,nuvtot
        vec(i) = rhs(i)
      enddo

      call gaussj (a, asav, nuvtot, nuvgauss, vec, 1, 1)

      return
      end

c-----------------------------------------------------------------------

      SUBROUTINE gaussj (a,asav,n,np,b,m,mp)

c     Solves A x = b using Gaussian Elimination
c     (from Numerical Recipes?). Returns solution in b.

      INTEGER m,mp,n,np,NMAX
      REAL a(np,np),b(np,mp), asav(np,np)
      PARAMETER (NMAX=500)
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      REAL big,dum,pivinv
      character*11 cline(900)

      do j=1,n
        do i=1,n
          asav(i,j) = a(i,j)
        enddo
      enddo

      do 11 j=1,n
        ipiv(j)=0
11    continue
      do 22 i=1,n
        big=0.
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (abs(a(j,k)).ge.big)then
                  big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
              else if (ipiv(k).gt.1) then
                pause 'singular matrix (i) in gaussj'
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
          do 15 l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.0.) then
          do iloop=1,2
            if (iloop.eq.1) iu = 6
            if (iloop.eq.2) iu = 20
            write (iu,*) 'a (singular ii) irow,icol=',irow,icol
            do io=1,n
              write (iu,'(a,i4,a,f15.8)')'i=',io,'  a(i,i)=',asav(io,io)
            enddo
            do iseg=1,(n-1)/100 + 1
             joa = (iseg-1)*100 + 1
             job = min (iseg*100, n)
             write (iu,'(/6x, 900i11)') (jo, jo=joa,job)
             do io=1,n
               do jo=joa,job
                 if (abs(asav(io,jo)).le.1.e-10) then
                   cline(jo) = '          Z'
                 else
                   write (cline(jo),'(f11.6)') asav(io,jo)
                 endif
                 if (io.eq.164) cline(jo)(2:2) = 'X'
                 if (io.eq.jo ) cline(jo)(2:2) = '*'
               enddo
               write (iu,'(i6, 900a11)') io, (cline(jo),jo=joa,job)
             enddo
             write (iu,'(6x, 900i11/)') (jo, jo=joa,job)
           enddo
          enddo
          pause 'singular (ii) matrix in gaussj'
        endif
        pivinv=1./a(icol,icol)
        a(icol,icol)=1.
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 17 l=1,m
          b(icol,l)=b(icol,l)*pivinv
17      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0.
            do 18 l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
            do 19 l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
19          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue
      return
      END
