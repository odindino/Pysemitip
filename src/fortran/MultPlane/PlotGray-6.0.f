c   ******************** PlotGray ********************** 
c                                                                       
c   for rendering and manipulating top views of images  
C
C   VERSION 6.0 - FEB/11
c                                                                       
      subroutine PlotGray(dat,ncdim,nrdim,ncol,nrow,ilog,iexpan,shade,
     &ibin,ixbin,iybin,px,py,pz,ipow,ycorrec,nfilt,ixfilt,iyfilt,irfilt,
     &disc,nscal,zsmin,zsmax,igraycut,grang,smin,smax,shad)
c      subroutine grey(ztop,shad,aux1,aux2,tmp1,tmp2,smin,smax,srang,          
c     &               iexpan,ish,jsh,ixbin,iybin,nx,ny,fn,auto,op)
c
      character*1 ztop(512,512),shad(512,512),tmp1(512,512),
     &tmp2(512,512),gcut,ans,buf(10000),op(2)
      real dat(ncdim,nrdim),aux1(262144),aux2(262144),plotx(256),
     &ploty(256),zsmin(10),zsmax(10)
      integer hist(256),ixfilt(10),iyfilt(10),irfilt(10),system
      logical auto
      pi=4.*atan(1.)
c
c   ilog=0 for regular plot, 1 for log plot
c   iexpan is expansion factor for plot
c   shade: 0=height->1=local height,2=curv,3=deriv,4=polar,5=d-filt,
c   6=variance,7=split
c   ibin is averaging size for curvature and derivative shading
c   ixbin,iybin averaging size for smooth background or local variance
c   px,py,pz is viewing vector for derivative shading
c   ipow is power for derivative shading
c   ycorrec is correction factor for polar and d-filt shading
c   nfilt,ixfilt,iyfilt,irfilt and filter parameters for polar shading
c   disc is discriminator for local variance
c   nscal,zsmin,zsmax are for split gray scales
c   igraycut=1 to accept optimal values, 0 to specify cuts, -1 to specify range
c   smin,smax are specified gray cuts
c   grang is specified gray scale range
c
      ixbin=1
      iybin=1
c      
c      if (.not.auto) write(6,*) 'xy expansion factor ?'
c      if (.not.auto) read(5,*) iexpan
c
c      shade=0                                                           
c      if (.not.auto)                                                    
c     &write(6,*) 'shade ? (0=height->1=local height,2=curv,3=deriv,',   
c     &'4=polar,5=d-filt,6=variance,7=split)'
c      if (.not.auto) read(5,*) shade                                    
      smin=1.e10                                                        
      smax=-1.e10                                                       
      dmin=1.e10                                                        
      dmax=-1.e10 
      if (shade.eq.7) go to 500
      if (shade.eq.6) go to 195                                         
      if (shade.eq.5) go to 170                                         
      if (shade.eq.4) go to 150                                         
      if (shade.eq.3) go to 130                                         
      if (shade.gt.1) go to 120                                         
      if (shade.gt.0) go to 210                                         
c                                                                       
c   height shading                                                      
c                                                                       
      nnrow=nrow                                                        
      nncol=ncol                                                        
      do 110 j=1,nnrow                                                  
         do 105 i=1,nncol                                               
            ij=i+(j-1)*ncol                                             
            aux1(ij)=dat(i,j)                                            
            aux2(ij)=dat(i,j)                                            
            if (aux1(ij).gt.dmax) dmax=aux1(ij)                         
            if (aux1(ij).lt.dmin) dmin=aux1(ij)                         
            if (aux2(ij).gt.smax) smax=aux2(ij)                         
            if (aux2(ij).lt.smin) smin=aux2(ij)                         
105      continue                                                       
110   continue                                                          
      go to 725                                                         
c                                                                       
c   curvature shading                                                   
c                                                                       
c120   write(6,*) 'average size (#pixels) ?'                             
c      read(5,*) ibin                                                    
120   if (ibin.le.1) ibin=1                                             
      ixbin=ibin                                                        
      iybin=ibin                                                        
      nnrow=nrow-(2*ibin+(ibin-1))                                      
      nncol=ncol-(2*ibin+(ibin-1))                                      
      ioff=ibin+((ibin-1)/2)                                            
      kmin=-((ibin-1)/2)                                                
      kmax=ibin/2                                                       
      do 125 j=1,nnrow                                                  
         do 123 i=1,nncol                                               
            ij=i+(j-1)*ncol                                             
            sum1=0.                                                     
            sum2=0.                                                     
            sum3=0.                                                     
            sum4=0.                                                     
            sum5=0.                                                     
            sum6=0.                                                     
            do 121 k=kmin,kmax                                          
               sum1=sum1+dat(i+ioff+k-ibin,j+ioff)                
               sum2=sum2+dat(i+ioff+k,j+ioff)                     
               sum3=sum3+dat(i+ioff+k+ibin,j+ioff)                
               sum4=sum4+dat(i+ioff,j+ioff+k-ibin)                
               sum5=sum5+dat(i+ioff,j+ioff+k)                     
               sum6=sum6+dat(i+ioff,j+ioff+k+ibin)                
121         continue                                                    
            aux1(ij)=dat(i+ioff,j+ioff)                            
            curv=-(sum1-2*sum2+sum3+sum4-2*sum5+sum6)/ibin              
            aux2(ij)=(2-shade)*aux1(ij)+(shade-1)*curv                    
            if (aux1(ij).gt.dmax) dmax=aux1(ij)                         
            if (aux1(ij).lt.dmin) dmin=aux1(ij)                         
            if (aux2(ij).gt.smax) smax=aux2(ij)                         
            if (aux2(ij).lt.smin) smin=aux2(ij)                         
123      continue                                                       
125   continue                                                          
      go to 725                                                         
c                                                                       
c   surface normal shading                                              
c                                                                       
c130   write(6,*) 'average size (#pixels) ?'                             
c      read(5,*) ibin                                                    
130   if (ibin.le.1) ibin=1                                             
      ixbin=ibin                                                        
      iybin=ibin                                                        
      nnrow=nrow-(ibin+(ibin-1))                                        
      nncol=ncol-(ibin+(ibin-1))                                        
      ioff=ibin                                                         
      kmin=0                                                            
      kmax=ibin-1                                                       
c      write(6,*) 'viewing vector (x,y,z) ?'                             
c      read(5,*) px,py,pz                                                
      pmag=sqrt(px**2+py**2+pz**2)                                      
      px=px/pmag                                                        
      py=py/pmag                                                        
      pz=pz/pmag                                                        
c      write(6,*) 'power (1->10) ?'                                      
c      read(5,*) ipow                                                    
      do 145 j=1,nnrow                                                  
         do 140 i=1,nncol                                               
            ij=i+(j-1)*ncol                                             
            sum1=0.                                                     
            sum2=0.                                                     
            sum3=0.                                                     
            sum4=0.                                                     
            do 135 k=kmin,kmax                                          
               sum1=sum1+dat(i+ioff+k-ibin,j+ioff)                
               sum2=sum2+dat(i+ioff+k,j+ioff)                     
               sum3=sum3+dat(i+ioff,j+ioff+k-ibin)                
               sum4=sum4+dat(i+ioff,j+ioff+k)                     
135         continue                                                    
            dx=-(sum2-sum1)/(ibin*ibin*xcalib)                          
            dy=-(sum4-sum3)/(ibin*ibin*xcalib)                          
            dz=1.                                                       
            dmag=sqrt(dx**2+dy**2+dz**2)                                
            dx=dx/dmag                                                  
            dy=dy/dmag                                                  
            dz=dz/dmag                                                  
            prod=dx*px+dy*py+dz*pz                                      
            aux1(ij)=dat(i+ioff,j+ioff)                            
            prod=((prod+1.)/2.)                                         
c           prod=abs(prod)                                              
c           aux2(ij)=(3-shade)*aux1(ij)+(shade-2)*prod**ipow              
            aux2(ij)=prod**ipow                                          
            if (aux1(ij).gt.dmax) dmax=aux1(ij)                         
            if (aux1(ij).lt.dmin) dmin=aux1(ij)                         
            if (aux2(ij).gt.smax) smax=aux2(ij)                         
            if (aux2(ij).lt.smin) smin=aux2(ij)                         
140      continue                                                       
145   continue                                                          
      go to 725                                                         
c                                                                       
c   2-d histogram of derivative values                                  
c                                                                       
150   do 154 j=1,512                                                    
         do 152 i=1,512                                                 
            aux1(i+(j-1)*512)=0.                                         
152      continue                                                       
154   continue                                                          
c      write(6,*) 'average size (#pixels) ?'                             
c      read(5,*) ibin                                                    
      if (ibin.le.1) ibin=1                                             
c      write(6,*) 'correction factor for y calibration ?'                
c      read(5,*) ycorrec                                                 
      ixbin=ibin                                                        
      iybin=ibin                                                        
      nnrow=nrow-(ibin+(ibin-1))                                        
      nncol=ncol-(ibin+(ibin-1))                                        
      ioff=ibin                                                         
      kmin=0                                                            
      kmax=ibin-1                                                       
      do 160 j=1,nnrow                                                  
         do 158 i=1,nncol                                               
            ij=i+(j-1)*ncol                                             
            sum1=0.                                                     
            sum2=0.                                                     
            sum3=0.                                                     
            sum4=0.                                                     
            do 155 k=kmin,kmax                                          
               sum1=sum1+dat(i+ioff+k-ibin,j+ioff)                
               sum2=sum2+dat(i+ioff+k,j+ioff)                     
               sum3=sum3+dat(i+ioff,j+ioff+k-ibin)                
               sum4=sum4+dat(i+ioff,j+ioff+k)                     
155         continue                                                    
            dx=-(sum2-sum1)/(ibin*ibin*xcalib)                          
            dy=-(sum4-sum3)/(ibin*ibin*xcalib*ycorrec)                  
            dz=1.                                                       
            dmag=sqrt(dx**2+dy**2+dz**2)                                
            dx=dx/dmag                                                  
            dy=dy/dmag                                                  
            dz=dz/dmag                                                  
c           dx=atan(dx)*(2./pi)                                         
c           dy=atan(dy)*(2./pi)                                         
            thetx=acos(dz)*(2./pi)*dx/sqrt(dx**2+dy**2)                 
            thety=acos(dz)*(2./pi)*dy/sqrt(dx**2+dy**2)                 
            ixp=nint((1.+iexpan*thetx)*256.)                            
            if (ixp.le.0) ixp=1                                         
            if (ixp.gt.512) ixp=512                                     
            iyp=nint((1.+iexpan*thety)*256.)                            
            if (iyp.le.0) iyp=1                                         
            if (iyp.gt.512) iyp=512                                     
            ixyp=ixp+(iyp-1)*512                                        
            aux1(ixyp)=aux1(ixyp)+1                                       
            if (aux1(ixyp).gt.smax) smax=aux1(ixyp)                       
158      continue                                                       
160   continue                                                          
      smin=0.                                                           
      nrow=512                                                          
      ncol=512                                                          
      nnrow=512                                                         
      nncol=512                                                         
      iexpan=1                                                          
      do 164 j=1,nrow                                                   
         do 162 i=1,ncol                                                
            ij=i+(j-1)*ncol                                             
            aux2(ij)=aux1(ij)                                             
162      continue                                                       
164   continue                                                          
      dmin=smin                                                         
      dmax=smax                                                         
      go to 725                                                         
c                                                                       
c   filter derivative image based on specified derivative values        
c                                                                       
c170   write(6,*) 'average size (#pixels) ?'                             
c      read(5,*) ibin                                                    
170   if (ibin.le.1) ibin=1                                             
c      write(6,*) 'correction factor for y calibration ?'                
c      read(5,*) ycorrec                                                 
      ixbin=ibin                                                        
      iybin=ibin                                                        
      nnrow=nrow-(ibin+(ibin-1))                                        
      nncol=ncol-(ibin+(ibin-1))                                        
      ioff=ibin                                                         
      kmin=0                                                            
      kmax=ibin-1                                                       
c      write(6,*) 'number of filters ?'                                  
c      read(5,*) nfilt                                                   
c      write(6,*) 'filter centers, radii ?'                              
c      do 172 k=1,nfilt                                                  
c         write(6,*) 'filter #',k                                        
c         write(6,*) 'ix,iy ?'                                           
c         read(5,*) ixfilt(k),iyfilt(k)                                  
c         write(6,*) 'radius ?'                                          
c         read(5,*) irfilt(k)                                            
c172   continue                                                          
      do 190 j=1,nnrow                                                  
         do 185 i=1,nncol                                               
            ij=i+(j-1)*ncol                                             
            aux1(ij)=dat(i+ioff,j+ioff)                            
            if (aux1(ij).gt.dmax) dmax=aux1(ij)                           
            if (aux1(ij).lt.dmin) dmin=aux1(ij)                           
            sum1=0.                                                     
            sum2=0.                                                     
            sum3=0.                                                     
            sum4=0.                                                     
            do 175 k=kmin,kmax                                          
               sum1=sum1+dat(i+ioff+k-ibin,j+ioff)                
               sum2=sum2+dat(i+ioff+k,j+ioff)                     
               sum3=sum3+dat(i+ioff,j+ioff+k-ibin)                
               sum4=sum4+dat(i+ioff,j+ioff+k)                     
175         continue                                                    
            dx=-(sum2-sum1)/(ibin*ibin*xcalib)                          
            dy=-(sum4-sum3)/(ibin*ibin*xcalib*ycorrec)                  
            dz=1.                                                       
            dmag=sqrt(dx**2+dy**2+dz**2)                                
            dx=dx/dmag                                                  
            dy=dy/dmag                                                  
            dz=dz/dmag                                                  
c           dx=atan(dx)*(2./pi)                                         
c           dy=atan(dy)*(2./pi)                                         
            thetx=acos(dz)*(2./pi)*dx/sqrt(dx**2+dy**2)                 
            thety=acos(dz)*(2./pi)*dy/sqrt(dx**2+dy**2)                 
            ixp=nint((1.+thetx)*256.)                                   
            iyp=nint((1.+thety)*256.)                                   
            aux2(ij)=0.                                                  
            do 180 k=1,nfilt                                            
               if (ixp.lt.(ixfilt(k)-irfilt(k))) go to 180              
               if (ixp.gt.(ixfilt(k)+irfilt(k))) go to 180              
               iytemp=sqrt(float(irfilt(k)**2-(ixp-ixfilt(k))**2))      
               if (iyp.lt.(iyfilt(k)-iytemp)) go to 180                 
               if (iyp.gt.(iyfilt(k)+iytemp)) go to 180                 
               aux2(ij)=1.                                               
180         continue                                                    
185      continue                                                       
190   continue                                                          
      smin=0.                                                           
      smax=1.                                                           
      go to 725                                                         
c                                                                       
c   compute local variance
c                                                                       
c195   write(6,*) 'region size x,y (#pixels) ?'                         
c      read(5,*) ixbin,iybin                                             
195   if (ixbin.lt.1) ixbin=1                                           
      if (iybin.lt.1) iybin=1
c      write(6,*) 'discriminator level (-1 for no discrimination) ?'
c      read(5,*) disc
      mcol=ixbin                                                        
      mrow=iybin                                                        
      nncol=ncol-mcol+1                                                 
      nnrow=nrow-mrow+1                                                 
      ioff=(mcol-1)/2                                                   
      joff=(mrow-1)/2
      sumv=0.
      do 205 j=1,nnrow                                                  
         do 203 i=1,nncol                                               
            ij=i+(j-1)*ncol                                             
            sum=0.                                                      
            do 197 jj=1,mrow                                            
               do 196 ii=1,mcol                                         
                  sum=sum+dat(i+ii-1,j+jj-2)                     
 196           continue                                                 
 197        continue                                                    
            aver=sum/float(mrow*mcol)
            sum=0.
            do 201 jj=1,mrow                                            
               do 200 ii=1,mcol                                         
                  sum=sum+(dat(i+ii-1,j+jj-2)-aver)**2
 200           continue                                                 
 201        continue                                                    
            var=sum/float(mrow*mcol)
            sumv=sumv+var
            aux1(ij)=dat(i+ioff,j+joff-1)                        
            if (disc.eq.-1.) then
               aux2(ij)=var
            else
               if (var.gt.disc) then
                  aux2(ij)=254.
               else
                  aux2(2)=0.
               end if
            end if
            if (aux1(ij).gt.dmax) dmax=aux1(ij)                           
            if (aux1(ij).lt.dmin) dmin=aux1(ij)                           
            if (aux2(ij).gt.smax) smax=aux2(ij)                           
            if (aux2(ij).lt.smin) smin=aux2(ij)                           
 203     continue                                                       
 205  continue 
      sumv=sumv/float(nnrow*nncol)
c      write(6,*) 'min,max,av variance =',smin,smax,sumv
      go to 725
c                                                                       
c   separate smooth background                                          
c                                                                       
c210   write(6,*) 'average size x,y (#pixels) ?'                         
c      read(5,*) ixbin,iybin                                             
210   if (ixbin.lt.1) ixbin=1                                           
      if (iybin.lt.1) iybin=1                                           
      mcol=ixbin                                                        
      mrow=iybin                                                        
      nncol=ncol-mcol+1                                                 
      nnrow=nrow-mrow+1                                                 
      ioff=(mcol-1)/2                                                   
      joff=(mrow-1)/2                                                   
      do 220 j=1,nnrow                                                  
         do 218 i=1,nncol                                               
            ij=i+(j-1)*ncol                                             
            sum=0.                                                      
            do 216 jj=1,mrow                                            
               do 214 ii=1,mcol                                         
                  sum=sum+dat(i+ii-1,j+jj-2)                     
214            continue                                                 
216         continue                                                    
            sum=sum/float(mrow*mcol)                                    
            aux1(ij)=dat(i+ioff,j+joff-1)                         
            aux2(ij)=aux1(ij)-shade*sum                                   
            if (aux1(ij).gt.dmax) dmax=aux1(ij)                           
            if (aux1(ij).lt.dmin) dmin=aux1(ij)                           
            if (aux2(ij).gt.smax) smax=aux2(ij)                           
            if (aux2(ij).lt.smin) smin=aux2(ij)                           
218      continue                                                       
220   continue
      go to 725
c
c   use split grey scale
c
 500  nnrow=nrow                                                        
      nncol=ncol                                                        
      do 520 j=1,nrow
         do 510 i=1,ncol
            ij=i+(j-1)*ncol
            aux1(ij)=dat(i,j)
            if (aux1(ij).gt.dmax) dmax=aux1(ij)                         
            if (aux1(ij).lt.dmin) dmin=aux1(ij)                         
 510     continue
 520  continue
c      write(6,*) 'min z, max z values =',dmin,dmax
c      write(6,*) 'number of grey-scales ?'
c      read(5,*) nscal
c      write(6,*) 'min z, max z positions for each scale'
c     &,' (in ascending order) ?'
c      do 525 is=1,nscal
c         read(5,*) zsmin(is),zsmax(is)
c 525  continue
      do 550 j=1,nrow
         do 540 i=1,ncol
            ij=i+(j-1)*ncol
            aux2(ij)=0.
            do 530 is=1,nscal
             if (dat(i,j).ge.zsmin(is).and.dat(i,j).lt.zsmax(is))
     &       aux2(ij)=1.+(dat(i,j)-zsmin(is))*254./(zsmax(is)-zsmin(is))
 530        continue
            if (aux2(ij).eq.0.and.dat(i,j).ge.zsmax(nscal))
     &       aux2(ij)=255.
            if (aux2(ij).gt.smax) smax=aux2(ij)                         
            if (aux2(ij).lt.smin) smin=aux2(ij)                         
 540     continue
 550  continue
c                                                                       
c   determine grey-scale range                                          
c                                                                       
725   dscal=254./(dmax-dmin)                                            
      do 728 j=1,nnrow                                                  
         do 727 i=1,nncol                                               
            ij=i+(j-1)*ncol                                             
            aux1(ij)=(aux1(ij)-dmin)*dscal                             
727      continue                                                       
728   continue                                                          
      sscal=254./(smax-smin)                                            
      srang=254./sscal                                                  
      do 730 i=1,256                                                    
         hist(i)=0                                                      
730   continue                                                          
      do 740 j=1,nnrow                                                  
         do 735 i=1,nncol                                               
            ij=i+(j-1)*ncol                                             
            aux2(ij)=(aux2(ij)-smin)*sscal                                
            if (aux2(ij).lt.0.) aux2(ij)=0.                               
            if (aux2(ij).gt.254.) aux2(ij)=254.
            if (ilog.eq.1) aux2(ij)=46.0*alog(aux2(ij)+1)-1.
            idat=aux2(ij)                                                
            hist(idat+1)=hist(idat+1)+1                                 
735      continue                                                       
740   continue                                                          
      do 745 i=1,256                                                    
         ploty(i)=hist(i)                                               
         plotx(i)=i-1.                                                  
745   continue                                                          
c                                                                       
c      if (auto) go to 750                                               
c
c      open(unit=12, file='graph.dat')
c      do 747 i=1,256,1
c          write(12,*) plotx(i),ploty(i)
c747   continue
c      close(unit=12)
c      isys=system('xgraph graph.dat &')
c                                                                       
750   isum1=0.                                                          
      do 755 i=1,256                                                    
         isum1=isum1+hist(i)                                            
755   continue                                                          
      isum1=isum1*0.02                                                  
      isum2=0.                                                          
      do 760 i=1,256                                                    
         isum2=isum2+hist(i)                                            
         if (isum2.gt.isum1) go to 770
760   continue
 770  smin=i-2
      isum2=0.                                                          
      do 780 i=256,1,-1                                                 
         isum2=isum2+hist(i)                                            
         if (isum2.gt.isum1) go to 790
780   continue
 790  smax=i
c      write(6,795) int(smin),int(smax),srang*(smax-smin)/254.
c 795  format(' grey-scale cut at ',2i5,'  (range = ',g12.3,
c     &' ang)')
c      ans='y'                                                           
c      if (.not.auto) write(6,*) 'ok (yes/no/specify) ?'
c      if (.not.auto) read(5,800) ans
c800   format(1a1)                                                       
c      if (ans.eq.'y') go to 820
      if (igraycut.eq.1) go to 820
c      if (ans.eq.'n') go to 815
      if (igraycut.eq.0) go to 818
c      write(6,*) 'desired grey-scale range ?'
c      read(5,*) grang
      srnew=(grang/srang)*254.
      sdif=srnew-(smax-smin)
      smin=smin-sdif/2.
      smax=smax+sdif/2.
      go to 818
c 815  write(6,*) 'window (min,max values) ?'                            
c      read(5,*) smin,smax                                               
818   write(6,*) ' grey-scale cut at ',int(smin),int(smax)
c                                                                       
820   sscal=254./(smax-smin)                                            
      srang=srang/sscal                                                 
c      write(6,*) 'grey-scale range =',srang                             
      do 840 j=1,nnrow                                                  
         do 835 i=1,nncol                                               
            ij=i+(j-1)*ncol                                             
            aux2(ij)=(aux2(ij)-smin)*sscal                                
            if (aux2(ij).lt.0.) aux2(ij)=0.                               
            if (aux2(ij).gt.254.) aux2(ij)=254.                           
835      continue                                                       
840   continue                                                          
c                                                                       
c   construct the image                                                 
c                                                                       
 850  ny=iexpan*nnrow
      if (ny.gt.512) ny=512
      nx=iexpan*nncol
      if (nx.gt.512) nx=512                                             
      ish=(512-nx)/2.                                                   
c      write(6,*) 'xshift =',ish                                         
      jsh=(512-ny)/2.                                                   
c      write(6,*) 'yshift =',jsh                                         
      do 880 j=512,1,-1                                                 
         jj=(j-jsh-1)/iexpan+1                                          
         do 870 i=512,1,-1                                              
            ii=(i-ish-1)/iexpan+1                                       
            if (iexpan.gt.1) go to 855
c
c   no expansion of data
c
            if (j.le.jsh.or.jj.gt.nnrow) go to 860                      
            if (i.le.ish.or.ii.gt.nncol) go to 860
            iijj=ii+(jj-1)*ncol                                         
            ztop(i,j)=char(nint(aux1(iijj)))                             
            shad(i,j)=char(nint(aux2(iijj)))
            go to 870
c
c   expansion of data using linear interpolation
c
 855        if (j.le.jsh.or.jj.ge.nnrow) go to 860                      
            if (i.le.ish.or.ii.ge.nncol) go to 860
            delx=float(i-ish-1)/iexpan-(ii-1)                           
            dely=float(j-jsh-1)/iexpan-(jj-1)                           
            iip1=ii+1                                                   
            jjp1=jj+1                                                   
            d00=aux1(ii+(jj-1)*ncol)                                    
            d01=aux1(ii+(jjp1-1)*ncol)                                  
            d10=aux1(iip1+(jj-1)*ncol)                                  
            d11=aux1(iip1+(jjp1-1)*ncol)                                
            ztop(i,j)=char(nint((1.-delx)*(1.-dely)*d00+                
     &      (1.-delx)*dely*d01+delx*(1.-dely)*d10+delx*dely*d11))       
            s00=aux2(ii+(jj-1)*ncol)                                    
            s01=aux2(ii+(jjp1-1)*ncol)                                  
            s10=aux2(iip1+(jj-1)*ncol)                                  
            s11=aux2(iip1+(jjp1-1)*ncol)                                
            shad(i,j)=char(nint((1.-delx)*(1.-dely)*s00+                
     &      (1.-delx)*dely*s01+delx*(1.-dely)*s10+delx*dely*s11))       
            go to 870                                                   
c                            
860         ztop(i,j)=char(255)                                           
            shad(i,j)=char(255)                                           
870      continue                                                       
880   continue                                                          
c
      if (shade.ne.6) return
c
c   transfer aux2 to data array
c
      do 920 j=1,nnrow                                                  
         do 910 i=1,nncol                                               
            ij=i+(j-1)*nncol                                             
            dat(i,j)=aux2(ij)
 910     continue                                                       
 920  continue 
      nrow=nnrow
      ncol=nncol
c
      return                                                            
      end                                                               
c                                                                       
c   output PGM file for IRFAN display
c                                                                       
      subroutine dispir(ifile,img)           
c                                                                       
      character*1 img(512,512),xvtop(15)
      integer system
      data xvtop/'P','5',' ','5','1','2',' ','5','1','2',' '
     &,'2','5','5',' '/
      xvtop(3)=char(10)
      xvtop(7)=char(10)
      xvtop(11)=char(10)
      xvtop(15)=char(10)
c                                                                       
      nrecl=15+512*512
      go to (101,102,103,104,105,106,107,108,109) ifile
101      open(unit=91,file='img1.PGM',access='direct',recl=nrecl)
         go to 110
102      open(unit=91,file='img2.PGM',access='direct',recl=nrecl)
         go to 110
103      open(unit=91,file='img3.PGM',access='direct',recl=nrecl)
         go to 110
104      open(unit=91,file='img4.PGM',access='direct',recl=nrecl)
         go to 110
105      open(unit=91,file='img5.PGM',access='direct',recl=nrecl)
         go to 110
106      open(unit=91,file='img6.PGM',access='direct',recl=nrecl)
         go to 110
107      open(unit=91,file='img7.PGM',access='direct',recl=nrecl)
         go to 110
108      open(unit=91,file='img8.PGM',access='direct',recl=nrecl)
         go to 110
109      open(unit=91,file='img9.PGM',access='direct',recl=nrecl)
110   continue
c
      write(91,rec=1) (xvtop(i),i=1,15),
     &((img(i,512-j+1),i=1,512),j=1,512)
      close(unit=91)
c
C      istat=system('i_view32 img1.PGM')
      return                                                            
      end
