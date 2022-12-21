classdef Square
   properties
      Width
      Height
   end
   properties (Dependent)
      Area
   end
   
   methods
       
      %constructor 
      function self = Square(w,h)
         self.Width = w;
         self.Height = h;
      end
      
      %defines how properties are defined
      function obj = set.Width(obj,value)
         if (value > 0)
            obj.Width = value;
         else
            error('Property value must be positive')
         end
      end
      
      function obj = set.Height(obj,value)
         if (value > 0)
            obj.Height = value;
         else
            error('Property value must be positive')
         end
      end
      
      %update every time the area
      function a = get.Area(obj)
         a = obj.Width * obj.Height;
      end
      
      
      %function
      function out = halfArea(self)
         out = (self.Area)/2;
      end
      
   end
end