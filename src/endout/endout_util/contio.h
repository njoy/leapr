
   subroutine contio(nin,nout,nscr,a,nb,nw)
   !--------------------------------------------------------------------
   ! Utility routine for ENDF coded and blocked binary tapes.
   ! Read, write, and/or convert one control record.
   ! Positive units are coded, negative ones are blocked binary.
   ! If any unit is zero, it is not used.
   !--------------------------------------------------------------------
   ! externals
   integer::nin,nout,nscr,nb,nw
   real(kr)::a(*)
   ! internals
   integer::inin,inout,inscr,i
   character(11)::field(2)

   !--input
   if (nin.lt.0) then
      inin=iabs(nin)
      read(inin) math,mfh,mth,nb,nw,(a(i),i=1,6)
   else if (nin.gt.0) then
      read(nin,'(6e11.0,i4,i2,i3,i5)') (a(i),i=1,6),math,mfh,mth,nsp
   endif

   !--output
   nb=0
   nw=6
   c1h=a(1)
   c2h=a(2)
   l1h=nint(a(3))
   l2h=nint(a(4))
   n1h=nint(a(5))
   n2h=nint(a(6))
   a(3)=l1h
   a(4)=l2h
   a(5)=n1h
   a(6)=n2h
   if (nout.eq.0.and.nscr.eq.0) return
   inout=iabs(nout)
   if (nout.lt.0) then
      write(inout) math,mfh,mth,nb,nw,(a(i),i=1,6)
      inout=0
   endif
   inscr=iabs(nscr)
   if (nscr.lt.0) then
      write(inscr) math,mfh,mth,nb,nw,(a(i),i=1,6)
      inscr=0
   endif
   if (nout.le.0.and.nscr.le.0) return

   !--special path to write blanks on *end* cards
   if (math.eq.-1) call atend(inout,inscr)
   if (math.eq.-1) return
   if (math.eq.0) call amend(inout,inscr)
   if (math.eq.0) return
   if (mfh.eq.0) call afend(inout,inscr)
   if (mfh.eq.0) return
   if (mth.eq.0) call asend(inout,inscr)
   if (mth.eq.0) return

   !--format the output
   call a11(c1h,field(1))
   call a11(c2h,field(2))
   if (nscr.gt.0) then
      write(nscr,'(2a11,4i11,i4,i2,i3,i5)') field(1),field(2),&
        l1h,l2h,n1h,n2h,math,mfh,mth,nsc
      nsc=nsc+1
      if (nsc.gt.99999) nsc=1
   endif
   if (nout.gt.0) then
      write(nout,'(2a11,4i11,i4,i2,i3,i5)') field(1),field(2),&
        l1h,l2h,n1h,n2h,math,mfh,mth,nsh
      nsh=nsh+1
      if (nsh.gt.99999) nsh=1
   endif
   return
   end subroutine contio



