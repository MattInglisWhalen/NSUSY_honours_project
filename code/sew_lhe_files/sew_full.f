      program decay
      implicit none
      integer i,j,k
      double precision P(0:4,30),wgt,aqcd,aqed,scale
      double precision P2(0:4,30),wgt2,aqcd2,aqed2,scale2
      double precision P3(0:4,30)
      integer ndecay,decayid
      integer ievent,ievent2
      integer lun1,lun2,lun3,anti
      integer nexternal, ic(7,30)
      integer nexternal2, ic2(7,30),ic3(7,30)
      integer colourp,colourd
      logical done,done2,firstpass
      character*140 buff2
      character*140 buff3
      character*140 filein
      character*140 decayprod
      character*140 fileout
      lun1=20
      lun2=21
      lun3=22
      write(*,*) "enter the particle to decay"
      read(*,*) decayid
      write(*,*) " input file name"
      read(*,*) filein
      write(*,*) " decay produc file"
      read(*,*) decayprod
      write(*,*) "output file name"
      read(*,*) fileout
      firstpass=.true.
      
      open(lun1,FILE=filein,STATUS='OLD')
      open(lun2,FILE=decayprod,STATUS='OLD')
      open(lun3,FILE=fileout,STATUS='REPLACE')
      done=.false.
      do while(.not. done) 
         anti=1
         call read_event(lun1,P,wgt,nexternal,ic,ievent,scale,aqcd,aqed,
     c        buff2,done,anti)
         if(firstpass) then
            call write_comments(lun3)
            firstpass = .false.
           
         endif
         ndecay=0
         do i=1,nexternal
            do k=0,4
               P3(k,i+ndecay) = P(k,i)
            enddo
            do k =1,7
               ic3(k,i+ndecay) = ic(k,i)
            enddo
               
            if(ic(1,i)*ic(1,i) .eq. decayid*decayid .and. 
     c                           ic(6,i) .eq. 1) then
              
               anti=ic(1,i)/abs( ic(1,i) )

               call read_event(lun2,P2,wgt2,nexternal2,ic2,ievent2,
     c              scale2,aqcd2,aqed2,buff3,done2,anti)
             
               if(.not. done2) then
                  ic3(6,i+ndecay)=2
                  
                  if(ic(4,i) .ne. 0) then
                     colourp = ic(4,i)
                  endif
                  if(ic(5,i) .ne. 0) then
                     colourd = ic(5,i)
                     
                  endif
                  
                  do j=2,nexternal2
                     call boostx(P2(0,j),P(0,i),P3(0,i+ndecay+j-1))
                     do k=1,7 
                        ic3(k,i+ndecay+j-1)= ic2(k,j)
                     enddo
                     if((ic2(4,j) .ne. 0) .and. (ic2(4,j) .eq. ic2(4,1)
     c                                                      ) ) then
                        ic3(4,i+ndecay+j-1) = colourp
                     endif
                     if((ic2(5,j) .ne.0).and. (ic2(5,j) .eq.ic2(5,1))
     c                    ) then
                        ic3(5,i+ndecay+j-1) = colourd
                        
                     endif
                        
                     ic3(2,i+ndecay+j-1) = ic2(2,j)+i+ndecay-1
                     ic3(3,i+ndecay+j-1) = ic2(2,j)+i+ndecay-1
                     
                  enddo
                  do j=i+1,nexternal
                     if(ic(2,j) > i) ic(2,j)=ic(2,j) + nexternal2-1
                     if(ic(3,j) >i) ic(3,j) = ic(3,j) + nexternal2-1
                  enddo
                  ndecay = ndecay +nexternal2-1
               else 
                  write(*,*) 'No more decay product'
               endif
            endif
         enddo
     
         call write_event(lun3,P3,wgt,nexternal+ndecay,ic3,ievent,scale,
     c        aqcd,aqed,buff2)
        
      enddo
      close(lun1)
      close(lun2)
      close(lun3)
      end
      
      subroutine boostx(p,q , pboost)
c
c This subroutine performs the Lorentz boost of a four-momentum.  The
c momentum p is assumed to be given in the rest frame of q.  pboost is
c the momentum p boosted to the frame in which q is given.  q must be a
c timelike momentum.
c
c input:
c       real    p(0:3)         : four-momentum p in the q rest  frame
c       real    q(0:3)         : four-momentum q in the boosted frame
c
c output:
c       real    pboost(0:3)    : four-momentum p in the boosted frame
c
      implicit none
      double precision p(0:4),q(0:4),pboost(0:4),pq,qq,m,lf

      double precision rZero
      parameter( rZero = 0.0d0 )


c
      qq = q(1)**2+q(2)**2+q(3)**2



      if ( qq.ne.rZero ) then
         pq = p(1)*q(1)+p(2)*q(2)+p(3)*q(3)
         m = sqrt(q(0)**2-qq)
         lf = ((q(0)-m)*pq/qq+p(0))/m
         pboost(0) = (p(0)*q(0)+pq)/m
         pboost(1) =  p(1)+q(1)*lf
         pboost(2) =  p(2)+q(2)*lf
         pboost(3) =  p(3)+q(3)*lf
         pboost(4) = p(4)
      else
         pboost(0) = p(0)
         pboost(1) = p(1)
         pboost(2) = p(2)
         pboost(3) = p(3)
         pboost(4) = p(4)
      endif
c
      return
      end
