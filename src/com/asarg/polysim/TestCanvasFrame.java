package com.asarg.polysim;




import javax.swing.*;
import javax.swing.border.EtchedBorder;
import java.awt.*;
import java.awt.event.*;
import java.net.URL;

public class TestCanvasFrame extends JFrame implements MouseWheelListener, MouseMotionListener, MouseListener,KeyListener, ComponentListener{

    TestCanvas tc;
    JToolBar stepControlToolBar = new JToolBar();
    ActionListener actionListener;
    TileEditorWindow tileEditorWindow = new TileEditorWindow(800,600);
    ControlButton next = new ControlButton("forward");
    ControlButton prev = new ControlButton("backward");
    ControlButton play = new ControlButton("play");
    ControlButton fastf = new ControlButton("fast-forward");
    ControlButton fastb = new ControlButton("fast-backward");
    IconButton optionButton = new IconButton();
    JMenuBar mainMenu = new JMenuBar();
    // Menu Bar Items
    JMenuItem newMenuItem = new JMenuItem("New Assembly");
    JMenuItem loadAssemblyMenuItem = new JMenuItem("Load");
    JMenuItem saveAssemblyMenuItem = new JMenuItem("Save");
    JMenuItem saveAsMenuItem = new JMenuItem("Save as...");
    JMenuItem closeMenuItem = new JMenuItem("Close");

    JMenuItem importTileSetMenuItem = new JMenuItem("Import Tile Set");
    JMenuItem tileSetEditorMenuItem = new JMenuItem("Tile Set Editor");

    JMenuItem undoMenuItem = new JMenuItem("Undo");
    JMenuItem redoMenuItem = new JMenuItem("Redo");

    JMenuItem seedCreatorMenuItem = new JMenuItem("Seed Creator");
    JMenuItem tileSystemOptionsMenuItem = new JMenuItem("Options");

    int width;
    int height;

    Point lastMouseXY = new Point(width,height);
    int dragCount = 0;

    Assembly assembly;
    Frontier frontier;
    PolyTile frontierTile;


    public TestCanvasFrame(int w, int h, final Assembly assembly)
    {
        this.assembly = assembly;
        width = w;
        height = h;
        frontier = this.assembly.calculateFrontier();
        setLayout(new BorderLayout());

        tc = new TestCanvas();
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        add(tc, BorderLayout.SOUTH);
        tc.setSize(width, height);

        addActionListeners();

        addToolBars();

        addMenuBars();

        pack();
        setVisible(true);
        drawGrid();
    }

    private void addToolBars(){
        stepControlToolBar.add(fastb);
        fastb.setPreferredSize(new Dimension(30, 25));
        stepControlToolBar.add(prev);
        prev.setPreferredSize(new Dimension(30, 25));
        stepControlToolBar.add(play);
        play.setPreferredSize(new Dimension(30, 25));
        stepControlToolBar.add(next);
        next.setPreferredSize(new Dimension(30, 25));
        stepControlToolBar.add(fastf);
        fastf.setPreferredSize(new Dimension(30, 25));

        optionButton.setText(String.valueOf('\uf013'));
        stepControlToolBar.add(optionButton);

        stepControlToolBar.setBorder(new EtchedBorder());
        stepControlToolBar.setLayout(new FlowLayout(FlowLayout.CENTER));
        add(stepControlToolBar);
    }

    private void addMenuBars(){
        JMenu fileMenu = new JMenu("File");
        newMenuItem = new JMenuItem("New Assembly");
        newMenuItem.addActionListener(actionListener);
        loadAssemblyMenuItem = new JMenuItem("Load");
        loadAssemblyMenuItem.addActionListener(actionListener);
        saveAssemblyMenuItem = new JMenuItem("Save");
        saveAssemblyMenuItem.addActionListener(actionListener);
        saveAsMenuItem = new JMenuItem("Save as...");
        saveAsMenuItem.addActionListener(actionListener);
        closeMenuItem = new JMenuItem("Close");
        closeMenuItem.addActionListener(actionListener);
        fileMenu.add(newMenuItem);
        fileMenu.add(loadAssemblyMenuItem);
        fileMenu.add(saveAssemblyMenuItem);
        fileMenu.add(saveAsMenuItem);
        fileMenu.addSeparator();
        fileMenu.add(closeMenuItem);

        JMenu toolsMenu = new JMenu("Tools");
        importTileSetMenuItem = new JMenuItem("Import Tile Set");
        importTileSetMenuItem.addActionListener(actionListener);
        tileSetEditorMenuItem = new JMenuItem("Tile Set Editor");
        tileSetEditorMenuItem.addActionListener(actionListener);
        toolsMenu.add(importTileSetMenuItem);
        toolsMenu.add(tileSetEditorMenuItem);

        JMenu editMenu = new JMenu("Edit");
        undoMenuItem = new JMenuItem("Undo");
        undoMenuItem.addActionListener(actionListener);
        redoMenuItem = new JMenuItem("Redo");
        redoMenuItem.addActionListener(actionListener);
        editMenu.add(undoMenuItem);
        editMenu.add(redoMenuItem);

        JMenu tileSystemMenu = new JMenu("Tile System");
        seedCreatorMenuItem = new JMenuItem("Seed Creator");
        seedCreatorMenuItem.addActionListener(actionListener);
        tileSystemOptionsMenuItem = new JMenuItem("Options");
        tileSystemOptionsMenuItem.addActionListener(actionListener);
        tileSystemMenu.add(seedCreatorMenuItem);
        tileSystemMenu.add(tileSystemOptionsMenuItem);

        mainMenu.add(fileMenu);
        mainMenu.add(editMenu);
        mainMenu.add(toolsMenu);
        mainMenu.add(tileSystemMenu);
        setJMenuBar(mainMenu);
    }

    private void addActionListeners(){
        addMouseWheelListener(this);
        addMouseListener(this);
        addMouseMotionListener(this);
        addKeyListener(this);
        addComponentListener(this);
        actionListener = new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                if(e.getSource().equals(next))
                {
                    if(!frontier.isEmpty()) {
                        assembly.attach();
                        frontier = assembly.calculateFrontier();
                        drawGrid();
                    }
                }else if(e.getSource().equals(prev))
                {
                    if(!assembly.getAttached().isEmpty()) {
                        assembly.detach();
                        tc.reset();
                        frontier = assembly.calculateFrontier();
                        drawGrid();
                    }
                }else if(e.getSource().equals(play)){
                    while(!frontier.isEmpty()){
                        assembly.attach();
                        tc.reset();
                        frontier = assembly.calculateFrontier();
                        drawGrid();
                        try {
                            Thread.sleep(1000);
                        }catch(InterruptedException ie) {
                        }
                    }
                }else if(e.getSource().equals(fastb)){
                    while(!assembly.getAttached().isEmpty()){
                        assembly.detach();
                    }
                    tc.reset();
                    frontier = assembly.calculateFrontier();
                    drawGrid();
                }else if(e.getSource().equals(fastf)){
                    while(!frontier.isEmpty()){
                        assembly.attach();
                        frontier = assembly.calculateFrontier();
                    }
                    tc.reset();
                    drawGrid();
                } else if (e.getSource().equals(newMenuItem)){
                    System.out.println("new assembly");
                } else if (e.getSource().equals(tileSetEditorMenuItem)){
                    tileEditorWindow.setVisible(true);
                }
            }
        };

        next.addActionListener(actionListener);
        prev.addActionListener(actionListener);
        play.addActionListener(actionListener);
        fastf.addActionListener(actionListener);
        fastb.addActionListener(actionListener);
    }


    public void drawGrid()
    {
        PlaceFrontierOnGrid();
        tc.drawGrid(assembly.Grid);
        repaint();
        RemoveFrontierFromGrid();

    }

    private PolyTile getFrontierPolyTile(){
        if(frontierTile==null){
            frontierTile = new PolyTile("");
            frontierTile.setColor("428bca");
            String glues[]= {null, null, null, null};
            frontierTile.addTile(0,0, glues);
            return frontierTile;
        }else return frontierTile;
    }

    private void PlaceFrontierOnGrid(){
        for (FrontierElement fe : frontier){
            assembly.Grid.put(fe.getLocation(), getFrontierPolyTile().getTile(0,0));
        }
    }

    private void RemoveFrontierFromGrid(){
        for (FrontierElement fe : frontier){
            assembly.Grid.remove(fe.getLocation());
        }
    }

    public void zoomInDraw()
    {
        int tileDiameter = tc.getTileDiameter();
        if(tileDiameter< width/2)
            tc.setTileDiameter((int)(tileDiameter*1.25));
        else return;

        tc.reset();
        drawGrid();
        repaint();
    }
    public void zoomOutDraw()
    {
        int tileDiameter = tc.getTileDiameter();
        if(tileDiameter > 10) {
            tc.setTileDiameter((int) (tileDiameter * .75));
        }
        else return;

        tc.reset();
        drawGrid();
        repaint();
    }

    @Override
    public void mouseWheelMoved(MouseWheelEvent e) {
        if(e.getWheelRotation() == 1)
        {
           zoomOutDraw();
        }
        if(e.getWheelRotation() == -1)
        {
           zoomInDraw();
        }
    }

    @Override
    public void mouseDragged(MouseEvent e) {
        tc.translateOffset(e.getX() - lastMouseXY.x, e.getY() - lastMouseXY.y);
        lastMouseXY=e.getPoint();
        tc.reset();
        drawGrid();
        repaint();
    }

    @Override
    public void mouseMoved(MouseEvent e) {

    }

    @Override
    public void mouseClicked(MouseEvent e) {

    }

    @Override
    public void mousePressed(MouseEvent e) {
        lastMouseXY = e.getPoint();

    }

    @Override
    public void mouseReleased(MouseEvent e) {
        dragCount = 0;

    }

    @Override
    public void mouseEntered(MouseEvent e) {

    }

    @Override
    public void mouseExited(MouseEvent e) {

    }

    @Override
    public void keyTyped(KeyEvent e) {

    }

    @Override
    public void keyPressed(KeyEvent e) {
        if(e.getKeyCode() == KeyEvent.VK_PAGE_UP)
        {
            zoomInDraw();
        }else if(e.getKeyCode() == KeyEvent.VK_PAGE_DOWN)
        {
            zoomOutDraw();
        }
        else if(e.getKeyCode() == KeyEvent.VK_RIGHT)
        {
            if(!frontier.isEmpty()) {
                assembly.attach();
                frontier = assembly.calculateFrontier();
                drawGrid();
            }
        }
        else if(e.getKeyCode() == KeyEvent.VK_LEFT)
        {
            assembly.detach();
            tc.reset();
            frontier = assembly.calculateFrontier();
            drawGrid();
        }
    }

    @Override
    public void keyReleased(KeyEvent e) {

    }


    @Override
    public void componentResized(ComponentEvent e) {
        remove(tc);
        add(tc);
    }

    @Override
    public void componentMoved(ComponentEvent e) {

    }

    @Override
    public void componentShown(ComponentEvent e) {

    }

    @Override
    public void componentHidden(ComponentEvent e) {

    }
}
