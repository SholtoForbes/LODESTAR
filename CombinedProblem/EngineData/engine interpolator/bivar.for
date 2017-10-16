
c       Modification: 
c       -Converted to .for from .f90
c       -Commented empty lines, and unused code
c       D.Preller 30/7/2011     
c
c       Original code with full comments can be found at:
c       http://orion.math.iastate.edu/burkardt/f_src/f_src.html
c     
      subroutine idbvip (md,ndp,xd,yd,zd,nip,xi,yi,zi,iwk,wk)
c
      implicit none
c
      integer ndp, nip
c
      real ap, bp, cp, dp
      integer iip, itipv, itpv, iwk(31*ndp + nip), jwipl, jwipt,
     $ jwit, jwit0, jwiwk, jwiwl, jwiwp, jwwpd, md, nl, nt, ntsc
      real p00, p01, p02, p03, p04, p05, p10, p11, p12, p13, p14, p20,
     $ p21, p22, p23, p30, p31, p32, p40, p41, p50, wk(8*ndp), x0,
     $ xd(ndp), xi(nip), xs1, xs2, y0, yd(ndp), yi(nip), ys1, ys2,
     $ zd(ndp), zi(nip)
c
      save /idlc/
      save /idpt/
c
      common /idlc/ itipv,xs1,xs2,ys1,ys2,ntsc(9)
      common /idpt/ itpv,x0,y0,ap,bp,cp,dp, 
     $              p00,p10,p20,p30,p40,p50,p01,p11,p21,p31,p41, 
     $              p02,p12,p22,p32,p03,p13,p23,p04,p14,p05
c
c  Error check.
c
      if ( md < 1 .or. md > 3 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IDBVIP - Fatal error!'
        write ( *, '(a)' ) '  Input parameter MD out of range.'
        stop
      end if
c 
      if ( ndp < 4 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IDBVIP - Fatal error!'
        write ( *, '(a)' ) '  Input parameter NDP out of range.'
        stop
      end if
c 
      if ( nip < 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IDBVIP - Fatal error!'
        write ( *, '(a)' ) '  Input parameter NIP out of range.'
        stop
      end if
c 
      if ( md == 1 ) then
        iwk(1) = ndp
      else
        if ( ndp /= iwk(1) ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'IDBVIP - Fatal error!'
          write ( *, '(a)' ) '  MD = 2 or 3 but NDP was 
     $                          changed since last call.'
          stop
        end if
      end if
c 
      if ( md <= 2 ) then
        iwk(3) = nip
      else
        if ( nip < iwk(3) ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'IDBVIP - Fatal error!'
          write ( *, '(a)' ) '  MD = 3 but NIP was changed 
     $                          since last call.'
          stop
        end if
      end if
c
c  Allocation of storage areas in the IWK array.
c
      jwipt = 16
      jwiwl = 6*ndp+1
      jwiwk = jwiwl
      jwipl = 24*ndp+1
      jwiwp = 30*ndp+1
      jwit0 = 31*ndp
      jwwpd = 5*ndp+1
c
c  Triangulate the XY plane.
c
      if ( md == 1 ) then
c
        call idtang ( ndp, xd, yd, nt, iwk(jwipt), nl, iwk(jwipl), 
     $                iwk(jwiwl), iwk(jwiwp), wk )
c
        iwk(5) = nt
        iwk(6) = nl
c
        if ( nt == 0 ) then
          return
        end if
c
      else
c
        nt = iwk(5)
        nl = iwk(6)
c
      end if
c
c  Locate all points at which interpolation is to be performed.
c
      if ( md <= 2 ) then
c
        itipv = 0
        jwit = jwit0
c
        do iip = 1, nip
c
          jwit = jwit+1
c
          call idlctn ( ndp, xd, yd, nt, iwk(jwipt), nl, iwk(jwipl), 
     $                  xi(iip), yi(iip), iwk(jwit), iwk(jwiwk), wk )
c
        end do
c
      end if
c
c  Estimate the partial derivatives at all data points.
c
      call idpdrv ( ndp, xd, yd, zd, nt, iwk(jwipt), wk, wk(jwwpd) )
c
c  Interpolate the ZI values.
c
      itpv = 0
      jwit = jwit0
c
      do iip = 1, nip
c
        jwit = jwit + 1
c
        call idptip (ndp, xd, yd, zd, nt, iwk(jwipt), nl, iwk(jwipl),wk,
     $               iwk(jwit), xi(iip), yi(iip), zi(iip) )
c
      end do
c 
      return
      end
      subroutine idgrid (xd,yd,nt,ipt,nl,ipl,nxi,nyi,xi,yi,ngp,igp)
c
c*******************************************************************************
c
c! IDGRID organizes grid points for surface fitting.
c
c
c  Discussion:
c
c    IDGRID sorts the points in ascending order of triangle numbers and 
      implicit none
c
      integer nl, nt, nxi, nyi, igp(nxi*nyi), il0, il0t3, ilp1, ilp1t3,
     $ insd, ip1, ip2, ip3, ipl(3*nl), ipt(3*nt), it0, it0t3, ixi,
     $ iximn, iximx, iyi, izi, jigp0, jigp1, jigp1i, jngp0, jngp1,
     $ l,  ngp(2*(nt+2*nl)), ngp0, ngp1, nl0, nt0, nxinyi
      real spdt, u1, u2, u3, v1, v2, v3, vpdt, x1, x2, x3,xd(*),xi(nxi),
     $ xii, ximn, ximx, xmn, xmx, y1, y2, y3, yd(*), yi(nyi), yii,
     $ yimn, yimx, ymn, ymx
c
c  Statement functions
c
      spdt(u1,v1,u2,v2,u3,v3) = (u1-u2)*(u3-u2)+(v1-v2)*(v3-v2)
c
      vpdt(u1,v1,u2,v2,u3,v3) = (u1-u3)*(v2-v3)-(v1-v3)*(u2-u3)
c
c  Preliminary processing
c
      nt0 = nt
      nl0 = nl
      nxinyi = nxi * nyi
      ximn = min ( xi(1), xi(nxi) )
      ximx = max ( xi(1), xi(nxi) )
      yimn = min ( yi(1), yi(nyi) )
      yimx = max ( yi(1), yi(nyi) )
c
c  Determine grid points inside the data area.
c
      jngp0 = 0
      jngp1 = 2*(nt0+2*nl0)+1
      jigp0 = 0
      jigp1 = nxinyi + 1
c
      do it0 = 1, nt0
c
        ngp0 = 0
        ngp1 = 0
        it0t3 = it0*3
        ip1 = ipt(it0t3-2)
        ip2 = ipt(it0t3-1)
        ip3 = ipt(it0t3)
        x1 = xd(ip1)
        y1 = yd(ip1)
        x2 = xd(ip2)
        y2 = yd(ip2)
        x3 = xd(ip3)
        y3 = yd(ip3)
        xmn = min ( x1, x2, x3 )
        xmx = max ( x1, x2, x3 )
        ymn = min ( y1, y2, y3 )
        ymx = max ( y1, y2, y3 )
        insd = 0
c
        do ixi = 1, nxi
c
          if ( xi(ixi) < xmn .or. xi(ixi) > xmx ) then
            if ( insd == 0 ) then
              cycle
            end if
            iximx = ixi-1
            go to 23
          end if
c
          if ( insd /= 1 ) then
            insd = 1
            iximn = ixi
          end if
c
        end do
c 
        if ( insd == 0 ) then
          go to 38
        end if
c
        iximx = nxi
c
23      continue
c
        do iyi = 1, nyi
c
          yii = yi(iyi)
c
          if ( yii < ymn .or. yii > ymx ) then
            go to 37
          end if
c
          do ixi = iximn, iximx
c
            xii = xi(ixi)
            l = 0
            if ( vpdt(x1,y1,x2,y2,xii,yii) ) 36,25,26
c
25          continue
c
            l = 1
26          continue
            if(vpdt(x2,y2,x3,y3,xii,yii))     36,27,28
27          continue
            l = 1
28          continue
            if(vpdt(x3,y3,x1,y1,xii,yii))     36,29,30
29          continue
            l = 1
30          continue
            izi = nxi*(iyi-1)+ixi
c
            if ( l == 1 ) go to 31
c
            ngp0 = ngp0+1
            jigp0 = jigp0+1
            igp(jigp0) = izi
            go to 36
c
31          continue
c 
            do jigp1i = jigp1,nxinyi
              if ( izi == igp(jigp1i) ) then
                go to 36
              end if
            end do
c 
            ngp1 = ngp1+1
            jigp1 = jigp1-1
            igp(jigp1) = izi
c
36          continue
c
          end do
c
37        continue
c
        end do
c
38      continue
c
        jngp0 = jngp0+1
        ngp(jngp0) = ngp0
        jngp1 = jngp1-1
        ngp(jngp1) = ngp1
c
      end do
c
c  Determine grid points outside the data area.
c  in semi-infinite rectangular area.
c
      do il0 = 1, nl0
c
        ngp0 = 0
        ngp1 = 0
        il0t3 = il0*3
        ip1 = ipl(il0t3-2)
        ip2 = ipl(il0t3-1)
        x1 = xd(ip1)
        y1 = yd(ip1)
        x2 = xd(ip2)
        y2 = yd(ip2)
c
        xmn = ximn
        xmx = ximx
        ymn = yimn
        ymx = yimx
c
        if ( y2 >= y1 ) then
          xmn = min ( x1, x2 )
        end if
c
        if ( y2 <= y1 ) then
          xmx = max ( x1, x2 )
        end if
c
        if ( x2 <= x1 ) then
          ymn = min ( y1, y2 )
        end if
c
        if ( x2 >= x1 ) then
          ymx = max ( y1, y2 )
        end if
c
        insd = 0
c
        do ixi = 1, nxi
c
          if ( xi(ixi) < xmn .or. xi(ixi) > xmx ) then
            if ( insd == 0 ) then
              go to 42
            end if
            iximx = ixi-1
            go to 43
          end if
c
          if ( insd /= 1 ) then
            insd = 1
            iximn = ixi
          end if
c
42        continue
c
        end do
c
        if ( insd == 0 ) go to 58
c
        iximx = nxi
c
43      continue
c
        do iyi = 1, nyi
c
          yii = yi(iyi)
          if(yii<ymn.or.yii>ymx)        go to 57
c
          do ixi = iximn,iximx
c
            xii = xi(ixi)
            l = 0
            if(vpdt(x1,y1,x2,y2,xii,yii))     46,45,56
   45           l = 1
   46           if(spdt(x2,y2,x1,y1,xii,yii))     56,47,48
   47           l = 1
   48           if(spdt(x1,y1,x2,y2,xii,yii))     56,49,50
   49           l = 1
   50           izi = nxi*(iyi-1)+ixi
            if(l==1)    go to 51
            ngp0 = ngp0+1
            jigp0 = jigp0+1
            igp(jigp0) = izi
            go to 56
c 
51          continue
c 
            do jigp1i = jigp1,nxinyi
              if(izi==igp(jigp1i))     go to 56
            end do
c 
53          continue
c
            ngp1 = ngp1+1
            jigp1 = jigp1-1
            igp(jigp1) = izi
c
56          continue
c
          end do
c
57        continue
c
        end do
c
58      continue
c
        jngp0 = jngp0+1
        ngp(jngp0) = ngp0
        jngp1 = jngp1-1
        ngp(jngp1) = ngp1
c
c  In semi-infinite triangular area.
c
60      continue
c
        ngp0 = 0
        ngp1 = 0
        ilp1 = mod(il0,nl0)+1
        ilp1t3 = ilp1*3
        ip3 = ipl(ilp1t3-1)
        x3 = xd(ip3)
        y3 = yd(ip3)
        xmn = ximn
        xmx = ximx
        ymn = yimn
        ymx = yimx
        if(y3>=y2.and.y2>=y1)   xmn = x2
        if(y3<=y2.and.y2<=y1)   xmx = x2
        if(x3<=x2.and.x2<=x1)   ymn = y2
        if(x3>=x2.and.x2>=x1)   ymx = y2
        insd = 0
c
        do ixi = 1, nxi
c
          if ( xi(ixi) < xmn .or. xi(ixi) > xmx ) then
            if(insd==0)   go to 62
            iximx = ixi-1
            go to 63
          end if
c
          if ( insd /= 1 ) then
            insd = 1
            iximn = ixi
          end if
c
62        continue
c
        end do
c
        if(insd==0)     go to 78
c
        iximx = nxi
c 
63      continue
c
        do iyi = 1, nyi
c 
          yii = yi(iyi)
          if(yii<ymn.or.yii>ymx)        go to 77
c 
          do ixi = iximn, iximx
c 
            xii = xi(ixi)
            l = 0
            if ( spdt(x1,y1,x2,y2,xii,yii) )     66,65,76
   65           l = 1
   66           if(spdt(x3,y3,x2,y2,xii,yii))     70,67,76
   67           l = 1
   70           izi = nxi*(iyi-1)+ixi
c 
            if ( l /= 1 ) then
              ngp0 = ngp0+1
              jigp0 = jigp0+1
              igp(jigp0) = izi
              go to 76
            end if
c 
            do jigp1i = jigp1, nxinyi
              if(izi==igp(jigp1i)) go to 76
            end do
c 
            ngp1 = ngp1+1
            jigp1 = jigp1-1
            igp(jigp1) = izi
c 
76          continue
c
          end do
c 
77        continue
c
        end do
c 
78      continue
c
        jngp0 = jngp0+1
        ngp(jngp0) = ngp0
        jngp1 = jngp1-1
        ngp(jngp1) = ngp1
c 
      end do
c 
      return
      end
      subroutine idlctn (ndp,xd,yd,nt,ipt,nl,ipl,xii,yii,iti,iwk,wk)
c
c*******************************************************************************
c
c! IDLCTN finds the triangle that contains a point.
c
c
c  Discusstion:
c
      implicit none
c
      integer ndp, nl, nt, i1, i2, i3, idp, idsc(9), il1, il1t3,
     $ il2, ip1, ip2, ip3, ipl(3*nl), ipt(3*nt), isc, it0, it0t3,
     $ iti, itipv, itsc, iwk(18*ndp), jiwk, jwk, nl0, nt0, ntl,
     $ ntsc, ntsci
      real spdt, u1, u2, u3, v1, v2, v3, vpdt, wk(8*ndp), x0, x1, x2,
     $ x3, xd(ndp), xii, xmn, xmx, xs1, xs2, y0, y1, y2, y3, yd(ndp),
     $ yii, ymn, ymx, ys1, ys2
c
      save /idlc/
c
      common /idlc/ itipv,xs1,xs2,ys1,ys2,ntsc(9)
c
c  Statement functions
c
      spdt(u1,v1,u2,v2,u3,v3) = (u1-u2)*(u3-u2)+(v1-v2)*(v3-v2)
c
      vpdt(u1,v1,u2,v2,u3,v3) = (u1-u3)*(v2-v3)-(v1-v3)*(u2-u3)
c
c  Preliminary processing
c
      nt0 = nt
      nl0 = nl
      ntl = nt0+nl0
      x0 = xii
      y0 = yii
c
c  Processing for a new set of data points
c
      if ( itipv/=0)      go to 30
c
c  Divide the x-y plane into nine rectangular sections.
c
      xmn = xd(1)
      xmx = xd(1)
      ymn = yd(1)
      ymx = yd(1)
      do idp = 2, ndp
        xmn = min ( xd(idp), xmn )
        xmx = max ( xd(idp), xmx )
        ymn = min ( yd(idp), ymn )
        ymx = max ( yd(idp), ymx )
      end do
c 
      xs1 = ( xmn + xmn + xmx ) / 3.0E+00
      xs2 = ( xmn + xmx + xmx ) / 3.0E+00
      ys1 = ( ymn + ymn + ymx ) / 3.0E+00
      ys2 = ( ymn + ymx + ymx ) / 3.0E+00
c
c  Determine and store in the iwk array, triangle numbers of
c  the triangles associated with each of the nine sections.
c
      ntsc(1:9) = 0
      idsc(1:9) = 0
c 
      it0t3 = 0
      jwk = 0
c
      do it0 = 1, nt0
c
        it0t3 = it0t3+3
        i1 = ipt(it0t3-2)
        i2 = ipt(it0t3-1)
        i3 = ipt(it0t3)
        xmn = min ( xd(i1), xd(i2), xd(i3) )
        xmx = max ( xd(i1), xd(i2), xd(i3) )
        ymn = min ( yd(i1), yd(i2), yd(i3) )
        ymx = max ( yd(i1), yd(i2), yd(i3) )
c
        if ( ymn <= ys1 ) then
          if(xmn<=xs1)                   idsc(1) = 1
          if(xmx>=xs1.and.xmn<=xs2)    idsc(2) = 1
          if(xmx>=xs2)                   idsc(3) = 1
        end if
c
        if ( ymx >= ys1 .and. ymn <= ys2 ) then
          if(xmn<=xs1)                   idsc(4) = 1
          if(xmx>=xs1.and.xmn<=xs2)    idsc(5) = 1
          if(xmx>=xs2)                   idsc(6) = 1
        end if
c
        if(ymx<ys2)                   go to 25
        if(xmn<=xs1)                   idsc(7) = 1
        if(xmx>=xs1.and.xmn<=xs2)    idsc(8) = 1
        if(xmx>=xs2)                   idsc(9) = 1
c
25      continue
c
        do isc = 1, 9
          if ( idsc(isc) /= 0 ) then
            jiwk = 9*ntsc(isc)+isc
            iwk(jiwk) = it0
            ntsc(isc) = ntsc(isc)+1
            idsc(isc) = 0
          end if
        end do
c
c  Store in the wk array the minimum and maximum of the X and
c  Y coordinate values for each of the triangle.
c
        jwk = jwk+4
        wk(jwk-3) = xmn
        wk(jwk-2) = xmx
        wk(jwk-1) = ymn
        wk(jwk)   = ymx
c
      end do
c
      go to 60
c
c  Check if in the same triangle as previous.
c
30     continue
c
      it0 = itipv
c
      if(it0>nt0)      go to 40
c
      it0t3 = it0*3
      ip1 = ipt(it0t3-2)
      x1 = xd(ip1)
      y1 = yd(ip1)
      ip2 = ipt(it0t3-1)
      x2 = xd(ip2)
      y2 = yd(ip2)
      if(vpdt(x1,y1,x2,y2,x0,y0) < 0.0E+00 )      go to 60
      ip3 = ipt(it0t3)
      x3 = xd(ip3)
      y3 = yd(ip3)
      if(vpdt(x2,y2,x3,y3,x0,y0) < 0.0E+00 )      go to 60
      if(vpdt(x3,y3,x1,y1,x0,y0) < 0.0E+00 )      go to 60
c
      iti = it0
      itipv = it0
c
      return
c
c  Check if on the same border line segment.
c
40     continue
c
      il1 = it0 / ntl
      il2 = it0-il1*ntl
      il1t3 = il1*3
      ip1 = ipl(il1t3-2)
      x1 = xd(ip1)
      y1 = yd(ip1)
      ip2 = ipl(il1t3-1)
      x2 = xd(ip2)
      y2 = yd(ip2)
      if(il2/=il1)      go to 50
      if(spdt(x1,y1,x2,y2,x0,y0) < 0.0E+00 )      go to 60
      if(spdt(x2,y2,x1,y1,x0,y0) < 0.0E+00 )      go to 60
      if(vpdt(x1,y1,x2,y2,x0,y0) > 0.0E+00 )      go to 60
c  
      iti = it0
      itipv = it0
c
      return
c
c  Check if between the same two border line segments.
c
50     continue
c
      if(spdt(x1,y1,x2,y2,x0,y0) > 0.0E+00 )      go to 60
c
      ip3 = ipl(3*il2-1)
      x3 = xd(ip3)
      y3 = yd(ip3)
c
      if ( spdt(x3,y3,x2,y2,x0,y0) <= 0.0E+00 )  then
        iti = it0
        itipv = it0
        return
      end if
c
c  Locate inside the data area.
c  Determine the section in which the point in question lies.
c
60     continue
c 
      isc = 1
c
      if ( x0 >= xs1 ) then
        isc = isc+1
      end if
c
      if ( x0 >= xs2 ) then
        isc = isc+1
      end if
c
      if ( y0 >= ys1 ) then
        isc = isc+3
      end if
c
      if ( y0 >= ys2 ) then
        isc = isc+3
      end if
c
c  Search through the triangles associated with the section.
c
      ntsci = ntsc(isc)
      if(ntsci<=0)      go to 70
      jiwk = -9+isc
c
      do itsc = 1, ntsci
c
        jiwk = jiwk+9
        it0 = iwk(jiwk)
        jwk = it0*4
        if(x0<wk(jwk-3))    go to 61
        if(x0>wk(jwk-2))    go to 61
        if(y0<wk(jwk-1))    go to 61
        if(y0>wk(jwk))      go to 61
        it0t3 = it0*3
        ip1 = ipt(it0t3-2)
        x1 = xd(ip1)
        y1 = yd(ip1)
        ip2 = ipt(it0t3-1)
        x2 = xd(ip2)
        y2 = yd(ip2)
        if(vpdt(x1,y1,x2,y2,x0,y0)<0.0E+00 )    go to 61
        ip3 = ipt(it0t3)
        x3 = xd(ip3)
        y3 = yd(ip3)
c
        if ( vpdt(x2,y2,x3,y3,x0,y0) >= 0.0E+00 ) then
c
          if ( vpdt(x3,y3,x1,y1,x0,y0) >= 0.0E+00 ) then
            iti = it0
            itipv = it0
            return
          end if
c
        end if
c
61      continue
c
      end do
c
c  Locate outside the data area.
c
70     continue
c
      do il1 = 1, nl0
c
        il1t3 = il1*3
        ip1 = ipl(il1t3-2)
        x1 = xd(ip1)
        y1 = yd(ip1)
        ip2 = ipl(il1t3-1)
        x2 = xd(ip2)
        y2 = yd(ip2)
        if(spdt(x2,y2,x1,y1,x0,y0)<0.0E+00 )    go to 72
        if(spdt(x1,y1,x2,y2,x0,y0)<0.0E+00 )    go to 71
        if(vpdt(x1,y1,x2,y2,x0,y0)>0.0E+00 )    go to 72
        il2 = il1
        go to 75
c
   71       continue
c
        il2 = mod(il1,nl0)+1
        ip3 = ipl(3*il2-1)
        x3 = xd(ip3)
        y3 = yd(ip3)
        if(spdt(x3,y3,x2,y2,x0,y0)<=0.0E+00 )    go to 75
c
   72       continue
c
      end do
c 
      it0 = 1
      iti = it0
      itipv = it0
c
      return
c
   75     continue
c
      it0 = il1*ntl+il2
      iti = it0
      itipv = it0
c
      return
      end
      subroutine idpdrv ( ndp, xd, yd, zd, nt, ipt, pd, wk )
c
c*******************************************************************************
c
c! IDPDRV estimates first and second partial derivatives at data points.
c
c
c  Parameters:
      implicit none
c
      integer ndp, nt
c
      real d12, d23, d31, dx1, dx2, dy1, dy2, dz1, dz2, dzx1, dzx2,
     $ dzy1, dzy2
      real, parameter :: epsln = 1.0E-06
      integer idp, ipt(3*nt), ipti(3), it, iv, jpd, jpd0, jpdmx,
     $ jpt, jpt0, nt0
      real pd(5*ndp), vpx, vpxx, vpxy, vpy, vpyx, vpyy, vpz, vpzmn,
     $ w1(3), w2(3), wi, wk(ndp), xd(ndp), xv(3), yd(ndp), yv(3),
     $ zd(ndp), zv(3), zxv(3), zyv(3)
c
c  Preliminary processing.
c
      nt0 = nt
c
c  Clear the PD array.
c
      jpdmx = 5*ndp
c 
      pd(1:jpdmx) = 0.0E+00
c 
      wk(1:ndp) = 0.0E+00
c
c  Estimate ZX and ZY.
c
      do it = 1, nt0
c 
        jpt0 = 3*(it-1)
c 
        do iv = 1, 3
          jpt = jpt0+iv
          idp = ipt(jpt)
          ipti(iv) = idp
          xv(iv) = xd(idp)
          yv(iv) = yd(idp)
          zv(iv) = zd(idp)
        end do
c 
        dx1 = xv(2)-xv(1)
        dy1 = yv(2)-yv(1)
        dz1 = zv(2)-zv(1)
        dx2 = xv(3)-xv(1)
        dy2 = yv(3)-yv(1)
        dz2 = zv(3)-zv(1)
        vpx = dy1*dz2-dz1*dy2
        vpy = dz1*dx2-dx1*dz2
        vpz = dx1*dy2-dy1*dx2
        vpzmn = abs(dx1*dx2+dy1*dy2)*epsln
c 
        if ( abs(vpz) > vpzmn ) then
c 
          d12 = sqrt((xv(2)-xv(1))**2+(yv(2)-yv(1))**2)
          d23 = sqrt((xv(3)-xv(2))**2+(yv(3)-yv(2))**2)
          d31 = sqrt((xv(1)-xv(3))**2+(yv(1)-yv(3))**2)
          w1(1) = 1.0E+00 / (d31*d12)
          w1(2) = 1.0E+00 / (d12*d23)
          w1(3) = 1.0E+00 / (d23*d31)
          w2(1) = vpz*w1(1)
          w2(2) = vpz*w1(2)
          w2(3) = vpz*w1(3)
c 
          do iv = 1, 3
            idp = ipti(iv)
            jpd0 = 5*(idp-1)
            wi = (w1(iv)**2)*w2(iv)
            pd(jpd0+1) = pd(jpd0+1)+vpx*wi
            pd(jpd0+2) = pd(jpd0+2)+vpy*wi
            wk(idp) = wk(idp)+vpz*wi
          end do
c 
        end if
c 
      end do
c 
      do idp = 1, ndp
        jpd0 = 5*(idp-1)
        pd(jpd0+1) = -pd(jpd0+1)/wk(idp)
        pd(jpd0+2) = -pd(jpd0+2)/wk(idp)
      end do
c
c  Estimate ZXX, ZXY, and ZYY.
c
      do it = 1, nt0
c 
        jpt0 = 3*(it-1)
c 
        do iv = 1, 3
          jpt = jpt0+iv
          idp = ipt(jpt)
          ipti(iv) = idp
          xv(iv) = xd(idp)
          yv(iv) = yd(idp)
          jpd0 = 5*(idp-1)
          zxv(iv) = pd(jpd0+1)
          zyv(iv) = pd(jpd0+2)
        end do
c 
        dx1 = xv(2)-xv(1)
        dy1 = yv(2)-yv(1)
        dzx1 = zxv(2)-zxv(1)
        dzy1 = zyv(2)-zyv(1)
        dx2 = xv(3)-xv(1)
        dy2 = yv(3)-yv(1)
        dzx2 = zxv(3)-zxv(1)
        dzy2 = zyv(3)-zyv(1)
        vpxx = dy1*dzx2-dzx1*dy2
        vpxy = dzx1*dx2-dx1*dzx2
        vpyx = dy1*dzy2-dzy1*dy2
        vpyy = dzy1*dx2-dx1*dzy2
        vpz = dx1*dy2-dy1*dx2
        vpzmn = abs(dx1*dx2+dy1*dy2)*epsln
c 
        if ( abs(vpz) > vpzmn ) then
c 
          d12 = sqrt((xv(2)-xv(1))**2+(yv(2)-yv(1))**2)
          d23 = sqrt((xv(3)-xv(2))**2+(yv(3)-yv(2))**2)
          d31 = sqrt((xv(1)-xv(3))**2+(yv(1)-yv(3))**2)
          w1(1) = 1.0E+00 /(d31*d12)
          w1(2) = 1.0E+00 /(d12*d23)
          w1(3) = 1.0E+00 /(d23*d31)
          w2(1) = vpz*w1(1)
          w2(2) = vpz*w1(2)
          w2(3) = vpz*w1(3)
c 
          do iv = 1, 3
            idp = ipti(iv)
            jpd0 = 5*(idp-1)
            wi = (w1(iv)**2)*w2(iv)
            pd(jpd0+3) = pd(jpd0+3)+vpxx*wi
            pd(jpd0+4) = pd(jpd0+4)+(vpxy+vpyx)*wi
            pd(jpd0+5) = pd(jpd0+5)+vpyy*wi
          end do
c 
        end if
c 
      end do
c 
      do idp = 1, ndp
        jpd0 = 5*(idp-1)
        pd(jpd0+3) = -pd(jpd0+3)/wk(idp)
        pd(jpd0+4) = -pd(jpd0+4)/(2.0*wk(idp))
        pd(jpd0+5) = -pd(jpd0+5)/wk(idp)
      end do
c 
      return
      end
      subroutine idptip(ndp,xd,yd,zd,nt,ipt,nl,ipl,pdd,iti,xii,yii,zii)
c
c*******************************************************************************
c
c! IDPTIP performs interpolation, determining a value of Z given X and Y.
c
c
c  Modified:
c
c    19 February 2001
c
c  Parameters:
      implicit none
c
      integer ndp, nl, nt
c
      real a, aa, ab, ac, act2, ad, adbc, ap, b, bb, bc, bdt2,
     $ bp, c, cc, cd, cp, csuv, d, dd, dlt, dp, dx, dy,
     $ g1, g2, h1, h2, h3
      integer i, idp, il1, il2, ipl(3*nl), ipt(3*nt), it0,
     $ iti, itpv, jipl, jipt, jpd, jpdd, kpd, ntl
      real lu, lv, p0, p00, p01, p02, p03, p04, p05, p1,
     $ p10, p11, p12, p13, p14, p2, p20, p21, p22, p23,
     $ p3, p30, p31, p32, p4, p40, p41, p5, p50, pd(15),
     $ pdd(5*ndp), thsv, thus, thuv, thxu, u, v, x(3), x0,
     $ xd(*), xii, y(3), y0, yd(*), yii, z(3), z0, zd(*), zii,
     $ zu(3), zuu(3), zuv(3), zv(3), zvv(3)
c
      save /idpt/
c
      common /idpt/ itpv,x0,y0,ap,bp,cp,dp, 
     $              p00,p10,p20,p30,p40,p50,p01,p11,p21,p31,p41, 
     $              p02,p12,p22,p32,p03,p13,p23,p04,p14,p05
c
c  Preliminary processing
c
      it0 = iti
      ntl = nt+nl
c
      if ( it0 > ntl ) then
        il1 = it0/ntl
        il2 = it0-il1*ntl
        if(il1==il2)      go to 40
        go to 60
      end if
c
c  Calculation of ZII by interpolation.
c  Check if the necessary coefficients have been calculated.
c
      if ( it0 == itpv )     go to 30
c
c  Load coordinate and partial derivative values at the vertexes.
c
      jipt = 3*(it0-1)
      jpd = 0
c 
      do i = 1, 3
c 
        jipt = jipt+1
        idp = ipt(jipt)
        x(i) = xd(idp)
        y(i) = yd(idp)
        z(i) = zd(idp)
        jpdd = 5*(idp-1)
c 
        do kpd = 1, 5
          jpd = jpd+1
          jpdd = jpdd+1
          pd(jpd) = pdd(jpdd)
        end do
c 
      end do
c
c  Determine the coefficients for the coordinate system
c  transformation from the XY system to the UV system and vice versa.
c
      x0 = x(1)
      y0 = y(1)
      a = x(2)-x0
      b = x(3)-x0
      c = y(2)-y0
      d = y(3)-y0
      ad = a*d
      bc = b*c
      dlt = ad-bc
      ap =  d/dlt
      bp = -b/dlt
      cp = -c/dlt
      dp =  a/dlt
c
c  Convert the partial derivatives at the vertexes of the
c  triangle for the UV coordinate system.
c
      aa = a*a
      act2 = 2.0E+00 *a*c
      cc = c*c
      ab = a*b
      adbc = ad+bc
      cd = c*d
      bb = b*b
      bdt2 = 2.0E+00 *b*d
      dd = d*d
c 
      do i = 1, 3
        jpd = 5*i
        zu(i) = a*pd(jpd-4)+c*pd(jpd-3)
        zv(i) = b*pd(jpd-4)+d*pd(jpd-3)
        zuu(i) = aa*pd(jpd-2)+act2*pd(jpd-1)+cc*pd(jpd)
        zuv(i) = ab*pd(jpd-2)+adbc*pd(jpd-1)+cd*pd(jpd)
        zvv(i) = bb*pd(jpd-2)+bdt2*pd(jpd-1)+dd*pd(jpd)
      end do
c
c  Calculate the coefficients of the polynomial.
c
      p00 = z(1)
      p10 = zu(1)
      p01 = zv(1)
      p20 = 0.5E+00 * zuu(1)
      p11 = zuv(1)
      p02 = 0.5E+00 * zvv(1)
      h1 = z(2)-p00-p10-p20
      h2 = zu(2)-p10-zuu(1)
      h3 = zuu(2)-zuu(1)
      p30 =  10.0E+00 * h1 - 4.0E+00 * h2 + 0.5E+00 * h3
      p40 = -15.0E+00 * h1 + 7.0E+00 * h2           - h3
      p50 =   6.0E+00 * h1 - 3.0E+00 * h2 + 0.5E+00 * h3
      h1 = z(3)-p00-p01-p02
      h2 = zv(3)-p01-zvv(1)
      h3 = zvv(3)-zvv(1)
      p03 =  10.0E+00 * h1 - 4.0E+00 * h2 + 0.5E+00 * h3
      p04 = -15.0E+00 * h1 + 7.0E+00 * h2    -h3
      p05 =   6.0E+00 * h1 - 3.0E+00 * h2 + 0.5E+00 * h3
      lu = sqrt(aa+cc)
      lv = sqrt(bb+dd)
      thxu = atan2(c,a)
      thuv = atan2(d,b)-thxu
      csuv = cos(thuv)
      p41 = 5.0E+00*lv*csuv/lu*p50
      p14 = 5.0E+00*lu*csuv/lv*p05
      h1 = zv(2)-p01-p11-p41
      h2 = zuv(2)-p11-4.0E+00 * p41
      p21 =  3.0E+00 * h1-h2
      p31 = -2.0E+00 * h1+h2
      h1 = zu(3)-p10-p11-p14
      h2 = zuv(3)-p11- 4.0E+00 * p14
      p12 =  3.0E+00 * h1-h2
      p13 = -2.0E+00 * h1+h2
      thus = atan2(d-c,b-a)-thxu
      thsv = thuv-thus
      aa =  sin(thsv)/lu
      bb = -cos(thsv)/lu
      cc =  sin(thus)/lv
      dd =  cos(thus)/lv
      ac = aa*cc
      ad = aa*dd
      bc = bb*cc
      g1 = aa * ac*(3.0E+00*bc+2.0E+00*ad)
      g2 = cc * ac*(3.0E+00*ad+2.0E+00*bc)
      h1 = -aa*aa*aa*(5.0E+00*aa*bb*p50+(4.0E+00*bc+ad)*p41) 
     $     -cc*cc*cc*(5.0E+00*cc*dd*p05+(4.0E+00*ad+bc)*p14)
      h2 = 0.5E+00 * zvv(2)-p02-p12
      h3 = 0.5E+00 * zuu(3)-p20-p21
      p22 = (g1*h2+g2*h3-h1)/(g1+g2)
      p32 = h2-p22
      p23 = h3-p22
      itpv = it0
c
c  Convert XII and YII to UV system.
c
30     continue
c
      dx = xii-x0
      dy = yii-y0
      u = ap*dx+bp*dy
      v = cp*dx+dp*dy
c
c  Evaluate the polynomial.
c
      p0 = p00+v*(p01+v*(p02+v*(p03+v*(p04+v*p05))))
      p1 = p10+v*(p11+v*(p12+v*(p13+v*p14)))
      p2 = p20+v*(p21+v*(p22+v*p23))
      p3 = p30+v*(p31+v*p32)
      p4 = p40+v*p41
      p5 = p50
      zii = p0+u*(p1+u*(p2+u*(p3+u*(p4+u*p5))))
      return
c
c  Calculation of ZII by extrapolation in the rectangle.
c  Check if the necessary coefficients have been calculated.
c
40     continue
c
      if(it0==itpv)     go to 50
c
c  Load coordinate and partial derivative values at the end
c  points of the border line segment.
c
      jipl = 3*(il1-1)
      jpd = 0
c 
      do i = 1, 2
c 
        jipl = jipl+1
        idp = ipl(jipl)
        x(i) = xd(idp)
        y(i) = yd(idp)
        z(i) = zd(idp)
        jpdd = 5*(idp-1)
c 
        do kpd = 1, 5
          jpd = jpd+1
          jpdd = jpdd+1
          pd(jpd) = pdd(jpdd)
        end do
c 
      end do
c
c  Determine the coefficients for the coordinate system
c  transformation from the XY system to the UV system
c  and vice versa.
c
      x0 = x(1)
      y0 = y(1)
      a = y(2)-y(1)
      b = x(2)-x(1)
      c = -b
      d = a
      ad = a*d
      bc = b*c
      dlt = ad-bc
      ap =  d/dlt
      bp = -b/dlt
      cp = -bp
      dp =  ap
c
c  Convert the partial derivatives at the end points of the
c  border line segment for the UV coordinate system.
c
      aa = a*a
      act2 = 2.0E+00 * a * c
      cc = c*c
      ab = a*b
      adbc = ad+bc
      cd = c*d
      bb = b*b
      bdt2 = 2.0E+00 * b * d
      dd = d*d
c
      do i = 1, 2
        jpd = 5*i
        zu(i) = a*pd(jpd-4)+c*pd(jpd-3)
        zv(i) = b*pd(jpd-4)+d*pd(jpd-3)
        zuu(i) = aa*pd(jpd-2)+act2*pd(jpd-1)+cc*pd(jpd)
        zuv(i) = ab*pd(jpd-2)+adbc*pd(jpd-1)+cd*pd(jpd)
        zvv(i) = bb*pd(jpd-2)+bdt2*pd(jpd-1)+dd*pd(jpd)
      end do
c
c  Calculate the coefficients of the polynomial.
c
      p00 = z(1)
      p10 = zu(1)
      p01 = zv(1)
      p20 = 0.5E+00 * zuu(1)
      p11 = zuv(1)
      p02 = 0.5E+00 * zvv(1)
c
      h1 = z(2)-p00-p01-p02
      h2 = zv(2)-p01-zvv(1)
      h3 = zvv(2)-zvv(1)
c
      p03 =  10.0E+00 * h1 - 4.0E+00*h2+0.5E+00*h3
      p04 = -15.0E+00 * h1 + 7.0E+00*h2    -h3
      p05 =   6.0E+00 * h1 - 3.0E+00*h2+0.5E+00*h3
c
      h1 = zu(2)-p10-p11
      h2 = zuv(2)-p11
c
      p12 =  3.0E+00*h1-h2
      p13 = -2.0E+00*h1+h2
      p21 = 0.0E+00
      p23 = -zuu(2)+zuu(1)
      p22 = -1.5E+00*p23
c
      itpv = it0
c
c  Convert XII and YII to UV system.
c
50     continue
c
      dx = xii-x0
      dy = yii-y0
      u = ap*dx+bp*dy
      v = cp*dx+dp*dy
c
c  Evaluate the polynomial.
c
      p0 = p00+v*(p01+v*(p02+v*(p03+v*(p04+v*p05))))
      p1 = p10+v*(p11+v*(p12+v*p13))
      p2 = p20+v*(p21+v*(p22+v*p23))
      zii = p0+u*(p1+u*p2)
c
      return
c
c  Calculation of ZII by extrapolation in the triangle.
c  Check if the necessary coefficients have been calculated.
c
60     continue
c
      if ( it0 /= itpv ) then
c
c  Load coordinate and partial derivative values at the vertex of the triangle.
c
        jipl = 3*il2-2
        idp = ipl(jipl)
        x0 = xd(idp)
        y0 = yd(idp)
        z0 = zd(idp)
        jpdd = 5*(idp-1)
c 
        do kpd = 1, 5
          jpdd = jpdd+1
          pd(kpd) = pdd(jpdd)
        end do
c
c  Calculate the coefficients of the polynomial.
c
        p00 = z0
        p10 = pd(1)
        p01 = pd(2)
        p20 = 0.5E+00*pd(3)
        p11 = pd(4)
        p02 = 0.5E+00*pd(5)
        itpv = it0
c
      end if
c
c  Convert XII and YII to UV system.
c
      u = xii-x0
      v = yii-y0
c
c  Evaluate the polynomial.
c
      p0 = p00+v*(p01+v*p02)
      p1 = p10+v*p11
      zii = p0+u*(p1+u*p20)
c 
      return
      end
      subroutine idsfft (md,ndp,xd,yd,zd,nxi,nyi,nzi,xi,yi,zi,iwk,wk)
c
c*******************************************************************************
c
c! IDSFFT fits a smooth surface Z(X,Y) given irregular (X,Y,Z) data.
c
c
c  Discussion:
c
      implicit none
c
      integer ndp, nxi, nyi, nzi
c
      real ap, bp, cp, dp
      integer il1, il2, iti, itpv, iwk(31*ndp + nxi*nyi), ixi, iyi,
     $ izi, jig0mn, jig0mx, jig1mn, jig1mx, jigp, jngp, jwigp,
     $ jwigp0, jwipl, jwipt, jwiwl, jwiwp, jwngp, jwngp0, jwwpd,
     $ md, ngp0, ngp1, nl, nngp, nt
      real p00, p01, p02, p03, p04, p05, p10, p11, p12, p13, p14, p20, 
     $ p21, p22, p23, p30, p31, p32, p40, p41, p50, wk(6*ndp), x0, 
     $ xd(ndp), xi(nxi), y0, yd(ndp), yi(nyi), zd(ndp), zi(nzi,nyi)
c
      save /idpt/
c
      common /idpt/ itpv,x0,y0,ap,bp,cp,dp, 
     $              p00,p10,p20,p30,p40,p50,p01,p11,p21,p31,p41, 
     $              p02,p12,p22,p32,p03,p13,p23,p04,p14,p05
c
c  Error check.
c
      if ( md < 1 .or. md > 3 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IDSFFT - Fatal error!'
        write(*,*)'  Input parameter MD out of range.'
        stop
      end if
c 
      if ( ndp < 4 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IDSFFT - Fatal error!'
        write ( *, '(a)' ) '  Input parameter NDP out of range.'
        stop
      end if
c 
      if ( nxi < 1 .or. nyi < 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IDSFFT - Fatal error!'
        write ( *, '(a)' ) '  Input parameter NXI or NYI out of range.'
        stop
      end if
c 
      if ( nxi > nzi ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IDSFFT - Fatal error!'
        write ( *, '(a)' ) '  Input parameter NZI is less than NXI.'
        stop
      end if
c 
      if ( md <= 1 ) then
c
        iwk(1) = ndp
c
      else
c 
        if ( ndp /= iwk(1) ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'IDSFFT - Fatal error!'
          write ( *, '(a)' ) '  MD = 2 or 3 but ndp was changed 
     $                          since last call.'
          stop
        end if
c
      end if
c
      if ( md <= 2 ) then
c
        iwk(3) = nxi
        iwk(4) = nyi
c    
      else
c 
        if ( nxi /= iwk(3) ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'IDSFFT - Fatal error!'
          write ( *, '(a)' ) 'MD = 3 but nxi was changed 
     $                        since last call.'
          stop
        end if
c 
        if ( nyi /= iwk(4) ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'IDSFFT - Fatal error!'
          write ( *, '(a)' ) '  MD = 3 but nyi was changed 
     $                          since last call.'
          stop
        end if
c
      end if
c
c  Allocation of storage areas in the IWK array.
c
      jwipt = 16
      jwiwl = 6*ndp+1
      jwngp0 = jwiwl-1
      jwipl = 24*ndp+1
      jwiwp = 30*ndp+1
      jwigp0 = 31*ndp
      jwwpd = 5*ndp+1
c
c  Triangulate the XY plane.
c
      if ( md == 1 ) then
c 
        call idtang ( ndp, xd, yd, nt, iwk(jwipt), nl, iwk(jwipl), 
     $                iwk(jwiwl), iwk(jwiwp), wk )
c
        iwk(5) = nt
        iwk(6) = nl
c 
        if ( nt == 0 ) then
          return
        end if
c 
      else
c
        nt = iwk(5)
        nl = iwk(6)
c
      end if
c
c  Sort output grid points in ascending order of the triangle
c  number and the border line segment number.
c 
      if ( md <= 2 ) then
c 
        call idgrid ( xd, yd, nt, iwk(jwipt), nl, iwk(jwipl), nxi, 
     $                nyi, xi, yi, iwk(jwngp0+1), iwk(jwigp0+1) )
c 
      end if
c
c  Estimate partial derivatives at all data points.
c
      call idpdrv ( ndp, xd, yd, zd, nt, iwk(jwipt), wk, wk(jwwpd) )
c
c  Interpolate the ZI values.
c
      itpv = 0
      jig0mx = 0
      jig1mn = nxi*nyi+1
      nngp = nt+2*nl
c 
      do jngp = 1, nngp
c
        iti = jngp
c
        if ( jngp > nt ) then
          il1 = (jngp-nt+1)/2
          il2 = (jngp-nt+2)/2
          if(il2>nl) then
            il2 = 1
          end if
          iti = il1*(nt+nl)+il2
        end if
c
        jwngp = jwngp0+jngp
        ngp0 = iwk(jwngp)
c
        if ( ngp0 /= 0 ) then
c
          jig0mn = jig0mx+1
          jig0mx = jig0mx+ngp0
c 
          do jigp = jig0mn, jig0mx
c
            jwigp = jwigp0+jigp
            izi = iwk(jwigp)
            iyi = (izi-1)/nxi+1
            ixi = izi-nxi*(iyi-1)
c
            call idptip ( ndp, xd, yd, zd, nt,iwk(jwipt),nl,iwk(jwipl), 
     $                    wk, iti, xi(ixi), yi(iyi), zi(ixi,iyi) )
c
          end do
c 
        end if
c
        jwngp = jwngp0+2*nngp+1-jngp
        ngp1 = iwk(jwngp)
c
        if ( ngp1 /= 0 ) then
c
          jig1mx = jig1mn-1
          jig1mn = jig1mn-ngp1
c 
          do jigp = jig1mn, jig1mx
c
            jwigp = jwigp0+jigp
            izi = iwk(jwigp)
            iyi = (izi-1)/nxi+1
            ixi = izi-nxi*(iyi-1)
c
            call idptip ( ndp, xd, yd, zd, nt,iwk(jwipt),nl,iwk(jwipl), 
     $                    wk, iti, xi(ixi), yi(iyi), zi(ixi,iyi) )
c
          end do
c
        end if
c 
      end do
c 
      return
      end
      subroutine idtang ( ndp, xd, yd, nt, ipt, nl, ipl, iwl, iwp, wk )
c
c*******************************************************************************
c
c! IDTANG performs triangulation.
c
c
c  Discussion:
c
      implicit none
c
      integer ndp
c
      real dsqf, dsqi, dsqmn
      real, parameter :: epsln = 1.0E-06
      integer idxchg, il, ilf, iliv, ilt3, ilvs, ip, ip1, ip1p1, ip2,
     $ ip3, ipl(6*ndp), ipl1, ipl2, iplj1, iplj2, ipmn1, ipmn2,
     $ ipt(6*ndp-15), ipt1, ipt2, ipt3, ipti, ipti1, ipti2, irep,
     $ it, it1t3, it2t3, itf(2), its, itt3, itt3r, iwl(18*ndp),
     $ iwp(ndp), ixvs, ixvspv, jl1, jl2, jlt3, jp, jp1, jp2,
     $ jpc, jpmn, jpmx, jwl, jwl1, jwl1mn, nl, nl0, nlf, nlfc,
     $ nlft2, nln, nlnt3, nlsh, nlsht3, nlt3
      integer, parameter :: nrep = 100
      integer nt, nt0, ntf, ntt3, ntt3p3
      real sp, spdt, u1, u2, u3, v1, v2, v3, vp, vpdt, wk(ndp), x1,
     $ x2, x3, xd(ndp), xdmp, y1, y2, y3, yd(ndp), ydmp
c
c  Statement functions
c
      dsqf(u1,v1,u2,v2) = (u2-u1)**2+(v2-v1)**2
      spdt(u1,v1,u2,v2,u3,v3) = (u2-u1)*(u3-u1)+(v2-v1)*(v3-v1)
      vpdt(u1,v1,u2,v2,u3,v3) = (v3-v1)*(u2-u1)-(u3-u1)*(v2-v1)
c
c  Preliminary processing
c
      if ( ndp < 4 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IDTANG - Fatal error!'
        write ( *, '(a)' ) '  Input parameter NDP out of range.'
        stop
      end if
c
c  Determine IPMN1 and IPMN2, the closest pair of data points.
c
      dsqmn = dsqf(xd(1),yd(1),xd(2),yd(2))
      ipmn1 = 1
      ipmn2 = 2
c 
      do ip1 = 1, ndp-1
c 
        x1 = xd(ip1)
        y1 = yd(ip1)
        ip1p1 = ip1+1
c 
        do ip2 = ip1p1, ndp
c 
          dsqi = dsqf(x1,y1,xd(ip2),yd(ip2))
c 
          if ( dsqi == 0.0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'IDTANG - Fatal error!'
            write ( *, '(a)' ) '  Two of the input data points 
     $                            are identical.'
            stop
          end if
c 
          if(dsqi<dsqmn) then
            dsqmn = dsqi
            ipmn1 = ip1
            ipmn2 = ip2
          end if
c 
        end do
c 
      end do
c
c  Compute the midpoint of the closest two data points.
c
      xdmp = (xd(ipmn1)+xd(ipmn2)) / 2.0E+00
      ydmp = (yd(ipmn1)+yd(ipmn2)) / 2.0E+00
c
c  Sort the other (NDP-2) data points in ascending order of
c  distance from the midpoint and store the sorted data point
c  numbers in the IWP array.
c
      jp1 = 2
c 
      do ip1 = 1, ndp
        if ( ip1 /= ipmn1 .and. ip1 /= ipmn2 ) then
          jp1 = jp1+1
          iwp(jp1) = ip1
          wk(jp1) = dsqf(xdmp,ydmp,xd(ip1),yd(ip1))
        end if
      end do
c 
      do jp1 = 3, ndp-1
c 
        dsqmn = wk(jp1)
        jpmn = jp1
c 
        do jp2 = jp1, ndp
          if(wk(jp2)<dsqmn) then
            dsqmn = wk(jp2)
            jpmn = jp2
          end if
        end do
c 
        its = iwp(jp1)
        iwp(jp1) = iwp(jpmn)
        iwp(jpmn) = its
        wk(jpmn) = wk(jp1)
c 
      end do
c
c  If necessary, modify the ordering in such a way that the
c  first three data points are not collinear.
c
      x1 = xd(ipmn1)
      y1 = yd(ipmn1)
      x2 = xd(ipmn2)
      y2 = yd(ipmn2)
c 
      do jp = 3, ndp
        ip = iwp(jp)
        sp = spdt(xd(ip),yd(ip),x1,y1,x2,y2)
        vp = vpdt(xd(ip),yd(ip),x1,y1,x2,y2)
        if ( abs(vp) > ( abs(sp) * epsln ) )   go to 37
      end do
c 
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'IDTANG - Fatal error!'
      write ( *, '(a)' ) '  All collinear data points.'
      stop
c 
   37     continue
c 
      if ( jp /= 3 ) then
c 
        jpmx = jp
c 
        do jpc = 4, jpmx
          jp = jpmx+4-jpc
          iwp(jp) = iwp(jp-1)
        end do
c 
        iwp(3) = ip
c 
      end if
c
c  Form the first triangle.  
c
c  Store point numbers of the vertexes of the triangle in the IPT array, 
c  store point numbers of the border line segments and the triangle number in
c  the IPL array.
c
      ip1 = ipmn1
      ip2 = ipmn2
      ip3 = iwp(3)
c 
      if ( vpdt(xd(ip1),yd(ip1),xd(ip2),yd(ip2),xd(ip3),yd(ip3)) 
     $           < 0.0E+00 ) then
        ip1 = ipmn2
        ip2 = ipmn1
      end if
c 
      nt0 = 1
      ntt3 = 3
      ipt(1) = ip1
      ipt(2) = ip2
      ipt(3) = ip3
      nl0 = 3
      nlt3 = 9
      ipl(1) = ip1
      ipl(2) = ip2
      ipl(3) = 1
      ipl(4) = ip2
      ipl(5) = ip3
      ipl(6) = 1
      ipl(7) = ip3
      ipl(8) = ip1
      ipl(9) = 1
c
c  Add the remaining data points, one by one.
c
      do jp1 = 4, ndp

        ip1 = iwp(jp1)
        x1 = xd(ip1)
        y1 = yd(ip1)
c
c  Determine the first invisible and visible border line segments, iliv and
c  ilvs.
c
        do il = 1, nl0
c
          ip2 = ipl(3*il-2)
          ip3 = ipl(3*il-1)
          x2 = xd(ip2)
          y2 = yd(ip2)
          x3 = xd(ip3)
          y3 = yd(ip3)
          sp = spdt(x1,y1,x2,y2,x3,y3)
          vp = vpdt(x1,y1,x2,y2,x3,y3)
c
          if ( il == 1 ) then
            ixvs = 0
            if(vp<=(abs(sp)*(-epsln)))   ixvs = 1
            iliv = 1
            ilvs = 1
            go to 53
          end if
c
          ixvspv = ixvs
c
          if ( vp <= (abs(sp)*(-epsln)) ) then
            ixvs = 1
            if(ixvspv==1)      go to 53
            ilvs = il
            if(iliv/=1)        go to 54
            go to 53
          end if
c
          ixvs = 0
c
          if ( ixvspv /= 0 ) then
            iliv = il
            if(ilvs/=1)        go to 54
          end if
c
53         continue
c
        end do
c
        if(iliv==1.and.ilvs==1)  ilvs = nl0
c
54       continue
c
        if(ilvs<iliv)  ilvs = ilvs+nl0
c
c  Shift (rotate) the IPL array to have the invisible border
c  line segments contained in the first part of the array.
c
55       continue
c 
        if ( iliv /= 1 ) then
c 
          nlsh = iliv-1
          nlsht3 = nlsh*3
c 
          do jl1 = 1,nlsht3
            jl2 = jl1+nlt3
            ipl(jl2) = ipl(jl1)
          end do
c 
          do jl1 = 1,nlt3
            jl2 = jl1+nlsht3
            ipl(jl1) = ipl(jl2)
          end do
c 
          ilvs = ilvs-nlsh
c 
        end if
c
c  Add triangles to the IPT array, 
c  update border line segments in the IPL array, 
c  set flags for the border line segments to be reexamined in the IWL array.
c
        jwl = 0
c 
        do il = ilvs, nl0
c 
          ilt3 = il*3
          ipl1 = ipl(ilt3-2)
          ipl2 = ipl(ilt3-1)
          it   = ipl(ilt3)
c
c  Add a triangle to the IPT array.
c
          nt0 = nt0+1
          ntt3 = ntt3+3
          ipt(ntt3-2) = ipl2
          ipt(ntt3-1) = ipl1
          ipt(ntt3)   = ip1
c
c  Update border line segments in the IPL array.
c
          if ( il == ilvs ) then
            ipl(ilt3-1) = ip1
            ipl(ilt3)   = nt0
          end if
c 
          if ( il == nl0 ) then
            nln = ilvs+1
            nlnt3 = nln*3
            ipl(nlnt3-2) = ip1
            ipl(nlnt3-1) = ipl(1)
            ipl(nlnt3)   = nt0
          end if
c
c  Determine the vertex that does not lie on the border
c  line segments.
c
          itt3 = it*3
          ipti = ipt(itt3-2)
c 
          if ( ipti == ipl1 .or. ipti == ipl2 ) then
            ipti = ipt(itt3-1)
            if ( ipti == ipl1 .or. ipti == ipl2 ) then
              ipti = ipt(itt3)
            end if
          end if
c
c  Check if the exchange is necessary.
c
          if ( idxchg(xd,yd,ip1,ipti,ipl1,ipl2) /= 0 ) then
c
c  Modify the IPT array.
c
            ipt(itt3-2) = ipti
            ipt(itt3-1) = ipl1
            ipt(itt3)   = ip1
            ipt(ntt3-1) = ipti
            if(il==ilvs)  ipl(ilt3) = it
            if(il==nl0.and.ipl(3)==it)      ipl(3) = nt0
c
c  Set flags in the IWL array.
c
            jwl = jwl+4
            iwl(jwl-3) = ipl1
            iwl(jwl-2) = ipti
            iwl(jwl-1) = ipti
            iwl(jwl)   = ipl2
c
          end if
c 
        end do
c 
        nl0 = nln
        nlt3 = nlnt3
        nlf = jwl/2
c 
        if ( nlf == 0 ) then
          go to 79
        end if
c
c  Improve triangulation.
c
        ntt3p3 = ntt3+3
c
        do irep = 1, nrep
c    
          do ilf = 1,nlf
c
            ipl1 = iwl(2*ilf-1)
            ipl2 = iwl(2*ilf)
c
c  Locate in the ipt array two triangles on both sides of
c  the flagged line segment.
c
            ntf = 0
c
            do itt3r = 3,ntt3,3
              itt3 = ntt3p3-itt3r
              ipt1 = ipt(itt3-2)
              ipt2 = ipt(itt3-1)
              ipt3 = ipt(itt3)
              if(ipl1/=ipt1.and.ipl1/=ipt2.and. ipl1/=ipt3)    go to 71
              if(ipl2/=ipt1.and.ipl2/=ipt2.and. ipl2/=ipt3)    go to 71
              ntf = ntf+1
              itf(ntf) = itt3/3
              if(ntf==2)     go to 72
71             continue
            end do
c
            if ( ntf < 2 )       go to 76
c
c  Determine the vertexes of the triangles that do not lie
c  on the line segment.
c
72           continue
c
            it1t3 = itf(1)*3
            ipti1 = ipt(it1t3-2)
            if(ipti1/=ipl1.and.ipti1/=ipl2)    go to 73
            ipti1 = ipt(it1t3-1)
c
            if ( ipti1 == ipl1 .or. ipti1 == ipl2 ) then
              ipti1 = ipt(it1t3)
            end if
c
73           continue
c
            it2t3 = itf(2)*3
            ipti2 = ipt(it2t3-2)
            if(ipti2/=ipl1.and.ipti2/=ipl2)    go to 74
            ipti2 = ipt(it2t3-1)
            if(ipti2/=ipl1.and.ipti2/=ipl2)    go to 74
            ipti2 = ipt(it2t3)
c
c  Check if the exchange is necessary.
c
74           continue
c
            if(idxchg(xd,yd,ipti1,ipti2,ipl1,ipl2)==0) then
              go to 76
            end if
c
c  Modify the IPT array.
c
            ipt(it1t3-2) = ipti1
            ipt(it1t3-1) = ipti2
            ipt(it1t3)   = ipl1
            ipt(it2t3-2) = ipti2
            ipt(it2t3-1) = ipti1
            ipt(it2t3)   = ipl2
c
c  Set new flags.
c
            jwl = jwl+8
            iwl(jwl-7) = ipl1
            iwl(jwl-6) = ipti1
            iwl(jwl-5) = ipti1
            iwl(jwl-4) = ipl2
            iwl(jwl-3) = ipl2
            iwl(jwl-2) = ipti2
            iwl(jwl-1) = ipti2
            iwl(jwl)   = ipl1
            do jlt3 = 3,nlt3,3
              iplj1 = ipl(jlt3-2)
              iplj2 = ipl(jlt3-1)
c
              if((iplj1==ipl1.and.iplj2==ipti2).or. 
     $           (iplj2==ipl1.and.iplj1==ipti2)) then
                                   ipl(jlt3) = itf(1)
              end if
c
              if((iplj1==ipl2.and.iplj2==ipti1).or. 
     $           (iplj2==ipl2.and.iplj1==ipti1)) then
                                  ipl(jlt3) = itf(2)
              end if
c
            end do
c
76           continue
c
          end do
c 
          nlfc = nlf
          nlf = jwl/2
c
c  Reset the IWL array for the next round.
c
          if ( nlf == nlfc ) go to 79
c
          jwl1mn = 2*nlfc+1
          nlft2 = nlf*2
c 
          do jwl1 = jwl1mn,nlft2
            jwl = jwl1+1-jwl1mn
            iwl(jwl) = iwl(jwl1)
          end do
c 
          nlf = jwl/2
c
        end do
c
79      continue
c
      end do
c
c  Rearrange the IPT array so that the vertexes of each triangle
c  are listed counter-clockwise.
c
      do itt3 = 3, ntt3, 3
c 
        ip1 = ipt(itt3-2)
        ip2 = ipt(itt3-1)
        ip3 = ipt(itt3)
c 
        if(vpdt(xd(ip1),yd(ip1),xd(ip2),yd(ip2),xd(ip3),yd(ip3))
     $              < 0.0E+00 ) then
          ipt(itt3-2) = ip2
          ipt(itt3-1) = ip1
        end if
c 
      end do
c 
      nt = nt0
      nl = nl0
c 
      return
      end
      function idxchg ( x, y, i1, i2, i3, i4 )
c
c*******************************************************************************
c
c! IDXCHG determines whether two triangles should be exchanged.
c
c
c  Discussion:
c
      implicit none
c
      real a1sq, a2sq, a3sq, a4sq, c1sq, c3sq
      real, parameter :: epsln = 1.0E-06
      integer i1, i2, i3, i4, idx, idxchg
      real s1sq, s2sq, s3sq, s4sq, u1, u2, u3, u4, x(*), x1,
     $ x2, x3, x4, y(*), y1, y2, y3, y4
c
c  Preliminary processing
c
      x1 = x(i1)
      y1 = y(i1)
      x2 = x(i2)
      y2 = y(i2)
      x3 = x(i3)
      y3 = y(i3)
      x4 = x(i4)
      y4 = y(i4)
c
      idx = 0
c 
      u3 = (y2-y3)*(x1-x3)-(x2-x3)*(y1-y3)
      u4 = (y1-y4)*(x2-x4)-(x1-x4)*(y2-y4)
c 
      if ( u3 * u4 > 0.0E+00 ) then
c 
        u1 = (y3-y1)*(x4-x1)-(x3-x1)*(y4-y1)
        u2 = (y4-y2)*(x3-x2)-(x4-x2)*(y3-y2)
c
        a1sq = (x1-x3)**2+(y1-y3)**2
        a4sq = (x4-x1)**2+(y4-y1)**2
        c1sq = (x3-x4)**2+(y3-y4)**2
        a2sq = (x2-x4)**2+(y2-y4)**2
        a3sq = (x3-x2)**2+(y3-y2)**2
        c3sq = (x2-x1)**2+(y2-y1)**2
c
        s1sq = u1*u1 / (c1sq*max(a1sq,a4sq))
        s2sq = u2*u2 / (c1sq*max(a2sq,a3sq))
        s3sq = u3*u3 / (c3sq*max(a3sq,a1sq))
        s4sq = u4*u4 / (c3sq*max(a4sq,a2sq))
c 
        if ( min ( s3sq, s4sq ) - min ( s1sq, s2sq ) > epsln ) then
          idx = 1
        end if
c 
      end if
c 
      idxchg = idx
c 
      return
      end
c      subroutine timestamp ( )
c
c*******************************************************************************
c
c! TIMESTAMP prints the current YMDHMS date as a time stamp.
c
c
c  Example:
c      implicit none
c
c      character ( len = 8 ) ampm
c      integer d
c      character ( len = 8 ) date
c      integer h
c      integer m
c      integer mm
c      character ( len = 9 ), parameter, dimension(12) :: month = (/ 
c     $  'January  ', 'February ', 'March    ', 'April    ', 
c     $  'May      ', 'June     ', 'July     ', 'August   ', 
c     $  'September', 'October  ', 'November ', 'December ' /)
c      integer n
c      integer s
c      character ( len = 10 )  time
c      integer values(8)
c      integer y
c      character ( len = 5 ) zone
c
c      call date_and_time ( date, time, zone, values )
c
c      y = values(1)
c      m = values(2)
c      d = values(3)
c      h = values(5)
c      n = values(6)
c      s = values(7)
c      mm = values(8)
c
c      if ( h < 12 ) then
c        ampm = 'AM'
c      else if ( h == 12 ) then
c        if ( n == 0 .and. s == 0 ) then
c          ampm = 'Noon'
c        else
c          ampm = 'PM'
c        end if
c      else
c        h = h - 12
c        if ( h < 12 ) then
c          ampm = 'PM'
c        else if ( h == 12 ) then
c          if ( n == 0 .and. s == 0 ) then
c            ampm = 'Midnight'
c          else
c            ampm = 'AM'
c          end if
c        end if
c      end if
c
c      write ( *, '(a,1x,i2,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' )
c     $  trim ( month(m) ), d, y, h, ':', n, ':', s, '.', mm, trim (ampm)
c
c      return
c      end
c