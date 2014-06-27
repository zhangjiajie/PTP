ArrayList Species_list = new ArrayList();
ArrayList Branch_list = new ArrayList();
ArrayList tb_list = new ArrayList();
boolean curr_lock = false;
int Xrange = 1000;
int Yrange = 1000;
int text_size = 15;
int radius = 10;

void setup(){
        size(Xrange, Yrange);
        smooth();
        textSize(text_size);
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
            float w  = float(string_coords[4]); 
            Branch br = new Branch(x1, y1, x2, y2, w);
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
        randomSeed(222);
        
}

void draw(){
        background(255,255,255);
        for (int i=0 ; i<Branch_list.size(); i++){
            Branch br = (Branch)Branch_list.get(i);
            br.draw();
        }
        
        for (int i=0; i<Species_list.size(); i++){
          Species sp = (Species)Species_list.get(i);
          sp.updatelock(false);
        }
        
        TextBox tb = new TextBox(mouseX, mouseY+radius+2);
        
        curr_lock = false;
        for (int i=0; i<Species_list.size(); i++){
          //println("sp" + i + "currlock before draw:" + curr_lock);
          Species sp = (Species)Species_list.get(i);
          sp.updatelock(curr_lock);
          curr_lock = sp.draw(mouseX, mouseY, tb);
          //println("sp" + i + "currlock after draw:" + curr_lock);
        }
        
        for(int i=0; i<tb_list.size(); i++){
            TextBox tbi = (TextBox)tb_list.get(i);
            tbi.draw();
        }
        
        tb.draw();
        stroke(1);
}


void keyPressed() {
    if (key == 'c' || key == 'C'){
        tb_list = new ArrayList();
        for (int i=0; i<Species_list.size(); i++){
            Species spe = (Species)Species_list.get(i);
            spe.reset_selected();
        }
    }
}

class Branch{
    float x1, y1, x2, y2, weight;
    
    Branch(float x1, float y1, float x2, float y2, float weight){
        this.x1 = x1;
        this.y1 = y1;
        this.x2 = x2;
        this.y2 = y2;
        this.weight = weight;
    }
    
    void rescale(float xmultiple, float ymultiple, float xshift, float yshift){
          this.x1 = (this.x1 + xshift) * xmultiple + 10;
          this.y1 = (this.y1 + yshift) * ymultiple + 10;
          this.x2 = (this.x2 + xshift) * xmultiple + 10;
          this.y2 = (this.y2 + yshift) * ymultiple + 10;
    }
    
    void draw(){
        strokeWeight(this.weight);
        line(this.x1, this.y1, this.x2, this.y2);
        strokeWeight(1);
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
    
    void draw_name(TextBox tb){
        if (this.mouseover_species){
            stroke(255, 0, 0);
            if (this.mouseover){
                tb.add(name);
            }
        }
    }
    
    void draw() {
        
        if (this.mouseover_species){
            strokeWeight(2);
            stroke(255, 0, 0);
        }else{
            strokeWeight(1);
            stroke(192,192,192);
        }
        ellipse((int)this.x, (int)this.y, this.radius, this.radius);
        stroke(0);
        strokeWeight(2);
    }
    
    void showName(TextBox tb){
        tb.add(name);
    }
    
}

class Species{
    ArrayList Taxon = new ArrayList();
    boolean mouseover_species = false;
    boolean locked = false;
    boolean selected = false;
    int r = (int)random(0, 255);
    int g = (int)random(0, 255);
    int b = (int)random(0, 255);

    Species(){
    }
    
    void reset_selected(){
        this.selected = false;
    }

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
    
    boolean draw(int mx, int my, TextBox tb){
        is_mouse_over(mx, my);
        strokeWeight(2);
        fill(this.r, this.g, this.b);
        if (this.selected){
            for (int i = 0; i< Taxon.size(); i++){
                Taxa taxa = (Taxa)Taxon.get(i);
                taxa.update(true, true);
                //taxa.draw_name(tb);
                taxa.draw();
            }
        }
        else if (mouseover_species){
            for (int i = 0; i< Taxon.size(); i++){
                Taxa taxa = (Taxa)Taxon.get(i);
                taxa.update(taxa.is_mouse_over(mx, my), true);
                taxa.draw_name(tb);
                taxa.draw();
                if (mousePressed){
                    taxa.showName(tb);
                    if (mouseButton == RIGHT){
                        this.selected = true;
                        boolean flagtb = true;
                        for (int j=0; j <tb_list.size(); j++){
                            TextBox tbj = (TextBox)tb_list.get(j);
                            if (tb == tbj){
                                flagtb = false;
                            }
                        }
                        if(flagtb){
                            tb_list.add(tb);
                        }
                    }
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

class TextBox{
    ArrayList taxa_names = new ArrayList();
    float x, y;
    TextBox(float x, float y){
        this.x = x;
        this.y = y;
    }
    
    void add(String name){
        boolean flag = true;
        for (int i = 0; i< this.taxa_names.size(); i++){
            String s = (String)this.taxa_names.get(i);
            if (s.equals(name)){
                flag = false;
                break;
            }
        }
        if(flag){
            taxa_names.add(name);
        }
    }
    
    float[] text_box_size(){
        float[] wh = new float[2];
        float maxwidth = 0.0;
        float height = 0.0;
        for(int i =0; i <this.taxa_names.size(); i++){
            String tname = (String)taxa_names.get(i);
            float sw = textWidth(tname);
            if(sw > maxwidth){
                maxwidth = sw;
            }
            height = height + text_size + 5;
        }
        wh[0] = maxwidth + 10;
        wh[1] = height + 10;
        return wh;
    }
    
    void draw(){
        if(this.taxa_names.size() > 0 ){
        float[] wh = this.text_box_size();
        float w = wh[0];
        float h = wh[1];
        //check right
        if((this.x + w) >= Xrange){
            float diff = this.x + w - Xrange;
            this.x = this.x - diff - 5;
        }
        if(this.x <=0){
            this.x = 5;
        }
        
        //check bottom
        if((this.y + h) >= Yrange){
            float diff = this.y + h - Yrange;
            this.y = this.y - diff -5;
        }
        if (this.y <= 0){
            this.y = 5;
        }
        
        //draw the box
        stroke(255, 0, 0);
        fill(255,250,205);
        rect(this.x, this.y, w, h);
        stroke(0);
        
        //draw text
        for(int i=0; i<this.taxa_names.size(); i++){
            String name = (String)this.taxa_names.get(i);
            float tx = this.x + 5;
            float ty = this.y + 10 + radius + text_size * i + 5 * i;  
            stroke(0);
            fill(0);
            text(name, tx, ty);
        }
        
        fill(255,255,255);
      }
    }
}
