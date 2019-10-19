function c = dot_t_v(a,na,b,nb,dim)
% a is a tensor array of physical dimension dim
% na is the number of elements in na
% b is a vector array of physical dimension dim
% nb is the number of elements in nb
% PRE: na == nb || na==1 || nb==1

   els = eldsize(1,dim);

   if na==1 && nb==1
      if dim == 1
         c = zeros(els,1);
         c(1) = a(1)*b(1);
         return
      elseif dim == 2
         c = zeros(els,1);
         c(1) = a(1)*b(1) + a(3)*b(2);
         c(2) = a(2)*b(1) + a(4)*b(2);
         return
      else
         c = zeros(els,1);
         c(1) = a(1)*b(1) + a(4)*b(2) + a(7)*b(3);
         c(2) = a(2)*b(1) + a(5)*b(2) + a(8)*b(3);
         c(3) = a(3)*b(1) + a(6)*b(2) + a(9)*b(3);
         return
      end
   end
   
   
   
   if nb==1
      if dim == 1
         c = zeros(els,na);
         for ii=1:na
            c(1,ii) = a(1,ii)*b(1);            
         end
         return
      elseif dim == 2
         c = zeros(els,na);
         for ii=1:na         
            c(1,ii) = a(1,ii)*b(1) + a(3,ii)*b(2);
            c(2,ii) = a(2,ii)*b(1) + a(4,ii)*b(2);
         end
         return
      else
         c = zeros(els,na);
         for ii=1:na
            c(1,ii) = a(1,ii)*b(1) + a(4,ii)*b(2) + a(7,ii)*b(3);
            c(2,ii) = a(2,ii)*b(1) + a(5,ii)*b(2) + a(8,ii)*b(3);
            c(3,ii) = a(3,ii)*b(1) + a(6,ii)*b(2) + a(9,ii)*b(3);
         end
         return
      end
   end
   
   
   
   if na==1
      if dim == 1
         c = zeros(els,nb);
         for ii=1:nb
            c(1,ii) = a(1)*b(1,ii);
         end
         return
      elseif dim == 2
         c = zeros(els,nb);
         for ii=1:nb         
            c(1,ii) = a(1)*b(1,ii) + a(3)*b(2,ii);
            c(2,ii) = a(2)*b(1,ii) + a(4)*b(2,ii);
         end
         return
      else
         c = zeros(els,nb);
         for ii=1:nb
            c(1,ii) = a(1)*b(1,ii) + a(4)*b(2,ii) + a(7)*b(3,ii);
            c(2,ii) = a(2)*b(1,ii) + a(5)*b(2,ii) + a(8)*b(3,ii);
            c(3,ii) = a(3)*b(1,ii) + a(6)*b(2,ii) + a(9)*b(3,ii);
         end
         return
      end
   end
   
   
   
   if na==nb
      if dim == 1
         c = zeros(els,nb);
         for ii=1:nb
            c(1,ii) = a(1,ii)*b(1,ii);
         end
         return
      elseif dim == 2
         c = zeros(els,nb);
         for ii=1:nb         
            c(1,ii) = a(1,ii)*b(1,ii) + a(3,ii)*b(2,ii);
            c(2,ii) = a(2,ii)*b(1,ii) + a(4,ii)*b(2,ii);
         end
         return
      else
         c = zeros(els,nb);
         for ii=1:nb
            c(1,ii) = a(1,ii)*b(1,ii) + a(4,ii)*b(2,ii) + a(7,ii)*b(3,ii);
            c(2,ii) = a(2,ii)*b(1,ii) + a(5,ii)*b(2,ii) + a(8,ii)*b(3,ii);
            c(3,ii) = a(3,ii)*b(1,ii) + a(6,ii)*b(2,ii) + a(9,ii)*b(3,ii);
         end
         return
      end
   end
   
   
   
   c = [];

   
   
end