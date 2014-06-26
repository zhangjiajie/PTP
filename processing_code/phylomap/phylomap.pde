ArrayList Species_list = new ArrayList();
ArrayList Branch_list = new ArrayList();
boolean curr_lock = false;
int Xrange = 800;
int Yrange = 800;

void setup(){
        size(Xrange, Yrange);
        smooth();
        String lines[] = loadStrings("xxx.taxa.txt");
        //println("there are " + lines.length + " lines");
        double minX = 999999;
        double minY = 999999;
        double maxX = -999999;
        double maxY = -999999;  
        for (int i=0; i < lines.length; i++) {
            Species s = new Species();
            String[] taxalist = split(lines[i], ',');
            for (int j=0; j < taxalist.length; j++){
              String[] xyname = split(taxalist[j], ' ');
              Taxa t = new Taxa(float(xyname[0]), float(xyname[1]), xyname[2]);
              s.add(t);
              if (float(xyname[0]) < minX){minX = float(xyname[0]);}
              if (float(xyname[0]) > maxX){maxX = float(xyname[0]);}
              if (float(xyname[1]) < minY){minY = float(xyname[1]);}
              if (float(xyname[1]) > maxY){maxY = float(xyname[1]);}
              
            }
            Species_list.add(s);
            //println(lines[i]);
        }
        
        String lines2[] = loadStrings("xxx.line.txt");
        for (int i=0;i<lines2.length;i++){
            String[] string_coords = split(lines2[i], ',');
            float x1 = float(string_coords[0]);
            float y1 = float(string_coords[1]);
            float x2 = float(string_coords[2]);
            float y2 = float(string_coords[3]);
            Branch br = new Branch(x1, y1, x2, y2);
            Branch_list.add(br);
            if (x1 < minX) {minX = x1;}
            if (x1 > maxX) {maxX = x1;}
            if (y1 < minY) {minY = y1;}
            if (y1 > maxY) {maxY = y1;}
            if (x2 < minX) {minX = x2;}
            if (x2 > maxX) {maxX = x2;}
            if (y2 < minY) {minY = y2;}
            if (y2 > maxY) {maxY = y2;}
        }
        
        double x_range = maxX - minX;
        double y_range = maxY - minY;
        double scale_range = x_range;
        if(y_range > x_range){scale_range = y_range;}
        double scale_factor = 0;
        if(scale_range == x_range){
          scale_factor = (Xrange-20)/scale_range ;
        }else{
          scale_factor = (Yrange-20)/scale_range ;
        }
        
        for (int i =0 ; i<Species_list.size(); i++){
          Species sp = (Species)Species_list.get(i);
          sp.rescale(scale_factor, scale_factor, -minX, -minY);
        }
        
        for (int i=0 ; i<Branch_list.size(); i++){
            Branch br = (Branch)Branch_list.get(i);
            br.rescale((float)scale_factor, (float)scale_factor, (float)-minX, (float)-minY);
        }
        
}

void draw(){
        background(220,250,250);
        for (int i=0; i<Species_list.size(); i++){
          Species sp = (Species)Species_list.get(i);
          sp.updatelock(false);
        }
        curr_lock = false;
        for (int i=0; i<Species_list.size(); i++){
          //println("sp" + i + "currlock before draw:" + curr_lock);
          Species sp = (Species)Species_list.get(i);
          sp.updatelock(curr_lock);
          curr_lock = sp.draw(mouseX, mouseY);
          //println("sp" + i + "currlock after draw:" + curr_lock);
        }
        
        stroke(1);
        
        for (int i=0 ; i<Branch_list.size(); i++){
            Branch br = (Branch)Branch_list.get(i);
            br.draw();
        }
        
}

class Branch{
    float x1, y1, x2, y2;
    
    Branch(float x1, float y1, float x2, float y2){
        this.x1 = x1;
        this.y1 = y1;
        this.x2 = x2;
        this.y2 = y2;
    }
    
    void rescale(float xmultiple, float ymultiple, float xshift, float yshift){
          this.x1 = (this.x1 + xshift) * xmultiple + 10;
          this.y1 = (this.y1 + yshift) * ymultiple + 10;
          this.x2 = (this.x2 + xshift) * xmultiple + 10;
          this.y2 = (this.y2 + yshift) * ymultiple + 10;
    }
    
    void draw(){
        line(this.x1, this.y1, this.x2, this.y2);
    }
    
}


class Taxa{
    double x, y;
    int radius = 10;
    String name = "";
    boolean mouseover = false;
    boolean mouseover_species = false;
    
    Taxa(double x, double y, String name){
        this.x = x;
        this.y = y;
        this.name = name;
    }

    void rescale(double xmultiple, double ymultiple, double xshift, double yshift){
          this.x = (this.x + xshift) * xmultiple + 10;
          this.y = (this.y + yshift) * ymultiple + 10;
    }

    boolean is_mouse_over(int mx, int my){
                int ix = (int)x;
                int iy = (int)y;
                
    return sqrt((ix-mx)*(ix-mx) + (iy-my)*(iy-my)) <= radius;
    }
    
    void update(boolean mouseover, boolean mouseover_species){
        this.mouseover = mouseover;
        this.mouseover_species = mouseover_species;
    }
    
    void draw() {
        if (this.mouseover_species){
            stroke(255, 0, 0);
            if (this.mouseover){
                fill(0, 102, 153);
                int tx = (int)x;
                int ty = (int)y - radius - 2;
                if (tx <= 0){tx = 30;}
                if (ty <= 0){ty = 30;} 
                text(name, tx, ty);
                fill(255, 255, 255);
            }
        }else{
            stroke(0);
        }
        ellipse((int)this.x, (int)this.y, this.radius, this.radius);
    }
    
    void showName(){
                fill(0, 102, 153);
                int tx = (int)x;
                int ty = (int)y - radius - 2;
                if (tx <= 0){tx = 30;}
                if (ty <= 0){ty = 30;} 
                text(name, tx, ty);
                fill(255, 255, 255);
    }
    
}

class Species{
    ArrayList Taxon = new ArrayList();
    boolean mouseover_species = false;
    boolean locked = false;

    Species(){}

    void add(Taxa t){
        Taxon.add(t);
    }
    
    void rescale(double xmultiple, double ymultiple, double xshift, double yshift){
          for (int i = 0; i< Taxon.size(); i++){
            Taxa taxa = (Taxa)Taxon.get(i);
            taxa.rescale(xmultiple, ymultiple, xshift, yshift);
          }
    }
    
    void updatelock(boolean l){
          this.locked = l;
    }

    boolean is_mouse_over(int mx, int my){
        if(locked){
            mouseover_species = false;
            return false;
        }
        for (int i = 0; i< Taxon.size(); i++){
            Taxa taxa = (Taxa)Taxon.get(i);
            if (taxa.is_mouse_over(mx, my)){
                mouseover_species = true;
                return true;
            }
        }
        mouseover_species = false;
        return false;
    }
    
    boolean draw(int mx, int my){
        is_mouse_over(mx, my);
        strokeWeight(2);
        if (mouseover_species){
            for (int i = 0; i< Taxon.size(); i++){
                Taxa taxa = (Taxa)Taxon.get(i);
                taxa.update(taxa.is_mouse_over(mx, my), true);
                taxa.draw();
                if (mousePressed){
                    taxa.showName();
                }
            }
        }else{
            for (int i = 0; i< Taxon.size(); i++){
                Taxa taxa = (Taxa)Taxon.get(i);
                taxa.update(false, false);
                taxa.draw();
            }
        }
        return mouseover_species;
    }
}
