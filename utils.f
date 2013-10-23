c subroutine prerr
      subroutine prerr(vza,vzb,vb)
      double complex vza
      double precision vzb,vb
      write (*,*) "error:",vza,vzb,vb
      end
