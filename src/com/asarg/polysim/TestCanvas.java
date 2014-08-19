package  com.asarg.polysim;

import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.util.*;
import java.util.List;

public class TestCanvas extends JPanel {


    private int tileDiameter = 30;
    private Point center;
    BufferedImage canvasBFI;
    Graphics2D cg2d;
    Dimension res = new Dimension();

    TestCanvas(int w, int h)
    {
        center = new Point(w/2, h/2);
        canvasBFI = new BufferedImage(w, h, BufferedImage.TYPE_INT_ARGB);
        cg2d = canvasBFI.createGraphics();
        cg2d.setColor(Color.black);
        res.setSize(w,h);
        cg2d.setClip(0,0,w,h);
    }
    @Override
    public void paintComponent(Graphics g) {
        g.setColor(Color.black);

        Graphics2D g2 = (Graphics2D) g;
        g.drawImage(canvasBFI, 0, 0, null);

    }

    public void drawPolyTile(PolyTile pt)
    {
        Drawer.clearGraphics(cg2d);
        List<Tile> lt = pt.tiles;
        if(lt.isEmpty())System.out.println("List empty in draw polytile");
        for(Tile t : lt)
        {

            cg2d.setColor(Color.black);
            cg2d.fillRect(t.getLocation().x*tileDiameter + center.x - tileDiameter/2 , -t.getLocation().y*tileDiameter + center.y - tileDiameter/2, tileDiameter,tileDiameter);
            if(t.getGlueE() != null && !t.getGlueE().isEmpty())
            {
                cg2d.setColor(Color.cyan);
                cg2d.fillRect(t.getLocation().x*tileDiameter + center.x - tileDiameter/2 + tileDiameter, -t.getLocation().y*tileDiameter + center.y - tileDiameter/2, tileDiameter,tileDiameter);
                System.out.println(t.getGlueE());
            }
            if(t.getGlueW() != null && !t.getGlueW().isEmpty())
            {
                cg2d.setColor(Color.cyan);
                cg2d.fillRect(t.getLocation().x*tileDiameter + center.x - tileDiameter/2 - tileDiameter, -t.getLocation().y*tileDiameter + center.y - tileDiameter/2, tileDiameter,tileDiameter);
            }
            if(t.getGlueS() != null && !t.getGlueS().isEmpty())
            {
                cg2d.setColor(Color.cyan);
                cg2d.fillRect(t.getLocation().x*tileDiameter + center.x - tileDiameter/2 , -t.getLocation().y*tileDiameter + center.y - tileDiameter/2 + tileDiameter, tileDiameter,tileDiameter);
            }
            if(t.getGlueN() != null && !t.getGlueN().isEmpty())
            {
                cg2d.setColor(Color.cyan);
                cg2d.fillRect(t.getLocation().x*tileDiameter + center.x - tileDiameter/2 , -t.getLocation().y*tileDiameter + center.y - tileDiameter/2  - tileDiameter, tileDiameter,tileDiameter);
            }
        }
         repaint();
    }

    public void drawGrid(HashMap<Point, Tile> hmpt)
    {

        Drawer.TileDrawer.drawTiles(cg2d, hmpt.entrySet(), tileDiameter, center);

       repaint();

    }
    @Override
    public Dimension getPreferredSize()
    {
        return new Dimension(res.width,res.height);
    }
}