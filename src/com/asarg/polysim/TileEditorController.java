package com.asarg.polysim;

import com.asarg.polysim.adapters.graphics.raster.Drawer;
import com.asarg.polysim.adapters.graphics.raster.EditorCanvas;
import javafx.application.Platform;
import javafx.beans.property.SimpleBooleanProperty;
import javafx.beans.property.SimpleIntegerProperty;
import javafx.beans.property.SimpleStringProperty;
import javafx.beans.value.ChangeListener;
import javafx.beans.value.ObservableValue;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.embed.swing.SwingNode;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.fxml.FXML;
import javafx.fxml.Initializable;
import javafx.scene.Node;
import javafx.scene.control.*;
import javafx.scene.control.cell.PropertyValueFactory;
import javafx.scene.control.cell.TextFieldTableCell;
import javafx.scene.input.KeyCode;
import javafx.scene.input.KeyCodeCombination;
import javafx.scene.input.KeyCombination;
import javafx.scene.input.KeyEvent;
import javafx.scene.layout.AnchorPane;
import javafx.scene.paint.Color;
import javafx.stage.FileChooser;
import javafx.stage.Stage;
import javafx.util.Callback;
import javafx.util.Pair;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;
import java.io.File;
import java.net.URL;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Map;
import java.util.ResourceBundle;

/**
 * Created by ericmartinez on 1/12/15.
 */
public class TileEditorController implements Initializable {

    TileConfiguration tc;
    EditorCanvas canvas;
    SimulationNode current;

    @FXML
    AnchorPane ptAnchorPane;
    @FXML
    SwingNode swingNodeCanvas;
    @FXML
    AnchorPane inspectorPane;
    @FXML
    TitledPane tileEditorPane;
    @FXML
    ListView<PolyTile> listview_polytiles;
    @FXML
    Accordion accordion;
    @FXML
    TextField field_tile_label;
    @FXML
    TextField field_north_glue;
    @FXML
    TextField field_east_glue;
    @FXML
    TextField field_south_glue;
    @FXML
    TextField field_west_glue;
    @FXML
    TextField field_glue_1;
    @FXML
    TextField field_glue_2;
    @FXML
    TextField field_strength;
    @FXML
    ColorPicker colorpicker_color;
    @FXML
    Button btn_delete_tile;
    @FXML
    Button btn_add_glue;
    @FXML
    Button btn_delete_row;
    @FXML
    Button btn_update_assembly;
    @FXML
    TableView<Glue> table_gluetable;
    @FXML
    MenuItem menu_export;
    @FXML
    MenuItem menu_delete;
    @FXML
    MenuItem context_menu_delete;
    @FXML
    MenuItem context_menu_new;
    @FXML
    MenuItem menu_clear;
    @FXML
    MenuItem menu_clear_glues;
    @FXML
    MenuItem menu_close;

    ObservableList<Glue> glueData;
    SimpleBooleanProperty updateable = new SimpleBooleanProperty(false);

    //Listeners
    ChangeListener<String> listener_label = new ChangeListener<String>() {
        @Override
        public void changed(ObservableValue<? extends String> observable, String oldValue, String newValue) {
            if(oldValue!=newValue) {
                Tile selected = canvas.getSelectedTileProperty().get();
                if (selected != null) {
                    canvas.getSelectedTileProperty().get().setLabel(newValue);
                    redrawPolyTile(selected);
                    updateable.set(true);

                    Platform.runLater(new Runnable() {
                        @Override
                        public void run() {
                            field_tile_label.requestFocus();
                            field_tile_label.positionCaret(field_tile_label.getText().length());
                        }
                    });
                }
            }
        }
    };

    ChangeListener<String> listener_north_glue = new ChangeListener<String>() {
        @Override
        public void changed(ObservableValue<? extends String> observable, String oldValue, String newValue) {
            Tile selected = canvas.getSelectedTileProperty().get();
            if(selected!=null) {
                canvas.getSelectedTileProperty().get().setGlueN(newValue);
                redrawPolyTile(selected);
                updateable.set(true);

                Platform.runLater(new Runnable() {
                    @Override
                    public void run() {
                        field_north_glue.requestFocus();
                        field_north_glue.positionCaret(field_north_glue.getText().length());
                    }
                });
            }
        }
    };

    ChangeListener<String> listener_south_glue = new ChangeListener<String>() {
        @Override
        public void changed(ObservableValue<? extends String> observable, String oldValue, String newValue) {
            Tile selected = canvas.getSelectedTileProperty().get();
            if(selected!=null) {
                canvas.getSelectedTileProperty().get().setGlueS(newValue);
                redrawPolyTile(selected);
                updateable.set(true);

                Platform.runLater(new Runnable() {
                    @Override
                    public void run() {
                        field_south_glue.requestFocus();
                        field_south_glue.positionCaret(field_south_glue.getText().length());
                    }
                });
            }
        }
    };

    ChangeListener<String> listener_east_glue = new ChangeListener<String>() {
        @Override
        public void changed(ObservableValue<? extends String> observable, String oldValue, String newValue) {
            Tile selected = canvas.getSelectedTileProperty().get();
            if(selected!=null) {
                canvas.getSelectedTileProperty().get().setGlueE(newValue);
                redrawPolyTile(selected);
                updateable.set(true);

                Platform.runLater(new Runnable() {
                    @Override
                    public void run() {
                        field_east_glue.requestFocus();
                        field_east_glue.positionCaret(field_east_glue.getText().length());
                    }
                });
            }
        }
    };

    ChangeListener<String> listener_west_glue = new ChangeListener<String>() {
        @Override
        public void changed(ObservableValue<? extends String> observable, String oldValue, String newValue) {
            Tile selected = canvas.getSelectedTileProperty().get();
            if(selected!=null) {
                canvas.getSelectedTileProperty().get().setGlueW(newValue);
                redrawPolyTile(selected);
                updateable.set(true);

                Platform.runLater(new Runnable() {
                    @Override
                    public void run() {
                        field_west_glue.requestFocus();
                        field_west_glue.positionCaret(field_west_glue.getText().length());
                    }
                });
            }
        }
    };

    public static class Glue{
        private SimpleStringProperty glueA;
        private SimpleStringProperty glueB;
        private SimpleStringProperty strength;

        private Glue(String g1, String g2, String str){
            glueA = new SimpleStringProperty(g1);
            glueB = new SimpleStringProperty(g2);
            strength = new SimpleStringProperty(str);
        }

        public String getGlueA(){ return glueA.get(); }

        public String getGlueB(){ return glueB.get(); }

        public String getStrength(){ return strength.get(); }

        public void setGlueA(String a){ glueA.set(a); }
        public void setGlueB(String b){ glueB.set(b); }
        public void setStrength(String str){
            try{
                int i = Integer.parseInt(str);
                if(i>=0) strength.set(str);
            }catch(NumberFormatException npe){

            }

        }
    }

    @Override
    public void initialize(URL fxmlFileLocation, ResourceBundle resources) {
        listview_polytiles.setCellFactory(new Callback<ListView<PolyTile>, ListCell<PolyTile>>() {
            @Override
            public ListCell<PolyTile> call(ListView<PolyTile> list) {
                final ListCell<PolyTile> ptCell = new PolyTileCell();
                return ptCell;
            }
        });
        if (this.tc != null) {
            listview_polytiles.setItems(tc.getTiletypes());
        }

        canvas = new EditorCanvas(800, 600);

        ptAnchorPane.widthProperty().addListener(new ChangeListener<Number>() {
            @Override
            public void changed(ObservableValue<? extends Number> observable, Number oldValue, Number newValue) {
                if(canvas.getHeight()>0)
                    canvas.resize(newValue.intValue(), canvas.getHeight());
                else
                    canvas.resize(newValue.intValue(), 600);
            }
        });

        ptAnchorPane.heightProperty().addListener(new ChangeListener<Number>() {
            @Override
            public void changed(ObservableValue<? extends Number> observable, Number oldValue, Number newValue) {
                if(canvas.getWidth()>0)
                    canvas.resize(canvas.getWidth(), newValue.intValue());
                else
                    canvas.resize(400, newValue.intValue());
            }
        });

        Platform.runLater(new Runnable() {
            @Override
            public void run(){
                canvas.resize((int) ptAnchorPane.getWidth(), (int) ptAnchorPane.getHeight());
            }
        });

        swingNodeCanvas.setContent(canvas);
        menu_delete.setDisable(true);
        context_menu_delete.setDisable(true);
        btn_update_assembly.managedProperty().bind(btn_update_assembly.visibleProperty());
        btn_update_assembly.visibleProperty().bind(updateable);
        disableTileData();
        addListeners();
        accordion.setExpandedPane(tileEditorPane);

        glueData = FXCollections.<Glue>observableArrayList();

        TableColumn col1 = table_gluetable.getColumns().get(0);
        TableColumn col2 = table_gluetable.getColumns().get(1);
        TableColumn col3 = table_gluetable.getColumns().get(2);

        col1.setCellValueFactory(new PropertyValueFactory<Glue, String>("glueA"));
        col2.setCellValueFactory(new PropertyValueFactory<Glue, String>("glueB"));
        col3.setCellValueFactory(new PropertyValueFactory<Glue, String>("strength"));
        col1.setCellFactory(TextFieldTableCell.forTableColumn());
        col2.setCellFactory(TextFieldTableCell.forTableColumn());
        col3.setCellFactory(TextFieldTableCell.forTableColumn());

        col1.setOnEditCommit(
                new EventHandler<TableColumn.CellEditEvent<Glue, String>>() {
                    @Override
                    public void handle(TableColumn.CellEditEvent<Glue, String> t) {
                        Glue g = ((Glue) t.getTableView().getItems().get(
                                t.getTablePosition().getRow())
                        );
                        tc.getGlueFunction().remove(new Pair<String, String>(g.getGlueA(), g.getGlueB()));
                        g.setGlueA(t.getNewValue());
                        Pair<String, String> p = new Pair<String, String>(g.getGlueA(), g.getGlueB());
                        tc.getGlueFunction().put(p, Integer.parseInt(g.getStrength()));
                        updateable.set(true);
                    }
                }
        );

        col2.setOnEditCommit(
                new EventHandler<TableColumn.CellEditEvent<Glue, String>>() {
                    @Override
                    public void handle(TableColumn.CellEditEvent<Glue, String> t) {
                        Glue g = ((Glue) t.getTableView().getItems().get(
                                t.getTablePosition().getRow())
                        );
                        tc.getGlueFunction().remove(new Pair<String, String>(g.getGlueA(), g.getGlueB()));
                        g.setGlueB(t.getNewValue());
                        Pair<String, String> p = new Pair<String, String>(g.getGlueA(), g.getGlueB());
                        tc.getGlueFunction().put(p, Integer.parseInt(g.getStrength()));
                        updateable.set(true);
                    }
                }
        );

        col3.setOnEditCommit(
                new EventHandler<TableColumn.CellEditEvent<Glue, String>>() {
                    @Override
                    public void handle(TableColumn.CellEditEvent<Glue, String> t) {
                        Glue g = ((Glue) t.getTableView().getItems().get(
                                t.getTablePosition().getRow())
                        );
                        tc.getGlueFunction().remove(new Pair<String, String>(g.getGlueA(), g.getGlueB()));
                        g.setStrength(t.getNewValue());
                        Pair<String, String> p = new Pair<String, String>(g.getGlueA(), g.getGlueB());
                        tc.getGlueFunction().put(p, Integer.parseInt(g.getStrength()));
                        updateable.set(true);
                    }
                }
        );

        table_gluetable.setItems(glueData);
        table_gluetable.setEditable(true);

        btn_delete_row.visibleProperty().set(false);

        table_gluetable.getSelectionModel().selectedItemProperty().addListener(new ChangeListener<Glue>() {
            @Override
            public void changed(ObservableValue<? extends Glue> observable, Glue oldValue, Glue newValue) {
                if(newValue!=null){
                    btn_delete_row.visibleProperty().set(true);
                }else{
                    btn_delete_row.visibleProperty().set(false);
                }
            }
        });
    }

    public void selectPolyTileIndex(int i){
        listview_polytiles.getSelectionModel().select(i);
    }

    public void setTileConfiguration(TileConfiguration tc){
        this.tc = tc;
        if(this.tc!=null){
            listview_polytiles.setItems(tc.getTiletypes());
            for (Map.Entry<Pair<String, String>, Integer> glF : tc.getGlueFunction().entrySet()) {
                String gLabelL = glF.getKey().getKey();
                String gLabelR = glF.getKey().getValue();
                int strength = glF.getValue();
                Glue g = new Glue(gLabelL, gLabelR, ""+strength);
                glueData.add(g);
            }
        }
    }

    public void setSimulationNode(SimulationNode current){
        this.current = current;
    }

    private void setTileData(Tile t){
        field_tile_label.setText(t.getLabel());
        field_north_glue.setText(t.getGlueN());
        field_east_glue.setText(t.getGlueE());
        field_south_glue.setText(t.getGlueS());
        field_west_glue.setText(t.getGlueW());
        colorpicker_color.setValue(Color.valueOf("0x"+ t.getParent().getColor()));
    }

    private void disableTileData(){
        field_tile_label.setDisable(true);
        field_north_glue.setDisable(true);
        field_east_glue.setDisable(true);
        field_west_glue.setDisable(true);
        field_south_glue.setDisable(true);
        colorpicker_color.setDisable(true);
        field_tile_label.setText(null);
        field_north_glue.setText(null);
        field_east_glue.setText(null);
        field_west_glue.setText(null);
        field_south_glue.setText(null);
        colorpicker_color.setValue(Color.WHITE);
        btn_delete_tile.setDisable(true);
    };

    private void enableTileData(){
        field_tile_label.setDisable(false);
        field_north_glue.setDisable(false);
        field_east_glue.setDisable(false);
        field_west_glue.setDisable(false);
        field_south_glue.setDisable(false);
        colorpicker_color.setDisable(false);
        btn_delete_tile.setDisable(false);
    };

    private void addTileListeners(){
        field_tile_label.textProperty().addListener(listener_label);
        field_north_glue.textProperty().addListener(listener_north_glue);
        field_south_glue.textProperty().addListener(listener_south_glue);
        field_east_glue.textProperty().addListener(listener_east_glue);
        field_west_glue.textProperty().addListener(listener_west_glue);
    }

    private void removeTileListeners(){
        field_tile_label.textProperty().removeListener(listener_label);
        field_north_glue.textProperty().removeListener(listener_north_glue);
        field_south_glue.textProperty().removeListener(listener_south_glue);
        field_east_glue.textProperty().removeListener(listener_east_glue);
        field_west_glue.textProperty().removeListener(listener_west_glue);
    }

    private void redrawPolyTile(Tile selected){
        canvas.drawPolyTile();
        int index = listview_polytiles.getSelectionModel().getSelectedIndex();
        PolyTile pt = canvas.getPolyTile();
        listview_polytiles.getItems().remove(index);
        listview_polytiles.getItems().add(index, pt);
        listview_polytiles.getSelectionModel().select(index);
        if(selected!=null) canvas.selectTile(selected);
    }

    public void deleteSelectedPolyTile(){
        int index = listview_polytiles.getSelectionModel().getSelectedIndex();
        if(index>=0) {listview_polytiles.getItems().remove(index); updateable.set(true);};
    }

    private void addListeners(){
        addTileListeners();

        listview_polytiles.getSelectionModel().selectedItemProperty().addListener(new ChangeListener<PolyTile>() {
            @Override
            public void changed(ObservableValue<? extends PolyTile> observable, PolyTile oldValue, PolyTile newValue) {
                if(newValue!=null) {
                    canvas.setPolyTile(newValue);
                    menu_delete.setDisable(false);
                    context_menu_delete.setDisable(false);
                    disableTileData();
                }else {
                    menu_delete.setDisable(true);
                    context_menu_delete.setDisable(true);
                    disableTileData();
                }
            }
        });

        listview_polytiles.setOnKeyReleased(new EventHandler<KeyEvent>() {
            final KeyCombination backspace = new KeyCodeCombination(KeyCode.BACK_SPACE);
            final KeyCombination delete = new KeyCodeCombination(KeyCode.DELETE);
            @Override
            public void handle(KeyEvent event) {
                if(backspace.match(event) || delete.match(event)){
                    deleteSelectedPolyTile();
                    updateable.set(true);
                }
            }
        });

        canvas.getSelectedTileProperty().addListener(new ChangeListener<Tile>() {
            @Override
            public void changed(ObservableValue<? extends Tile> observable, Tile oldValue, final Tile newValue) {
                if(newValue==null){
                    Platform.runLater(new Runnable() {
                        @Override
                        public void run() {
                            disableTileData();
                        }
                    });
                }else {
                    Platform.runLater(new Runnable() {
                        @Override
                        public void run() {
                            removeTileListeners();
                            setTileData(newValue);
                            enableTileData();
                            addTileListeners();
                        }
                    });

                }
            }
        });

        colorpicker_color.setOnAction(new EventHandler<ActionEvent>() {
            @Override
            public void handle(ActionEvent event) {
                PolyTile pt = listview_polytiles.getSelectionModel().getSelectedItem();
                if(pt!=null) {
                    Color color  = colorpicker_color.getValue();
                    pt.setColor(String.format( "%02X%02X%02X",
                            (int)( color.getRed() * 255 ),
                            (int)( color.getGreen() * 255 ),
                            (int)( color.getBlue() * 255 ) ));
                    redrawPolyTile(canvas.getSelectedTileProperty().get());
                    updateable.set(true);

                    Platform.runLater(new Runnable() {
                        @Override
                        public void run() {
                            colorpicker_color.requestFocus();
                        }
                    });
                }
            }
        });

        btn_delete_tile.setOnAction(new EventHandler<ActionEvent>() {
            @Override
            public void handle(ActionEvent event) {
                Tile selectedTile = canvas.getSelectedTileProperty().get();
                canvas.getPolyTile().removeTile(selectedTile.getLocation().x, selectedTile.getLocation().y);
                redrawPolyTile(null);
                updateable.set(true);
            }
        });

        btn_add_glue.setOnAction(new EventHandler<ActionEvent>() {
            @Override
            public void handle(ActionEvent event) {
                String glue1 = field_glue_1.getText().trim();
                String glue2 = field_glue_2.getText().trim();
                if(!glue1.isEmpty() && !glue2.isEmpty()){
                    Glue g = new Glue(glue1, glue2, field_strength.getText().trim());
                    table_gluetable.getItems().add(g);
                    Pair<String, String> p = new Pair<String, String>(g.getGlueA(), g.getGlueB());
                    tc.getGlueFunction().put(p, Integer.parseInt(g.getStrength()));
                    field_glue_1.clear();
                    field_glue_2.clear();
                    field_strength.clear();
                    updateable.set(true);
                }

            }
        });

        btn_delete_row.setOnAction(new EventHandler<ActionEvent>() {
            @Override
            public void handle(ActionEvent event) {
                int index = table_gluetable.getSelectionModel().getSelectedIndex();
                Glue g = table_gluetable.getItems().get(index);
                table_gluetable.getItems().remove(index);
                Pair<String, String> p = new Pair<String, String>(g.getGlueA(), g.getGlueB());
                tc.getGlueFunction().remove(p);
                updateable.set(true);
            }
        });


        menu_delete.setOnAction(new EventHandler<ActionEvent>() {
            @Override
            public void handle(ActionEvent event) {
                deleteSelectedPolyTile();
            }
        });

        context_menu_new.setOnAction(new EventHandler<ActionEvent>() {
            @Override
            public void handle(ActionEvent event) {
                PolyTile toAdd = new PolyTile();
                listview_polytiles.getItems().add(0,toAdd);
                listview_polytiles.getSelectionModel().select(0);
            }
        });

        context_menu_delete.setOnAction(new EventHandler<ActionEvent>() {
            @Override
            public void handle(ActionEvent event) {
                deleteSelectedPolyTile();
            }
        });

        menu_clear.setOnAction(new EventHandler<ActionEvent>() {
            @Override
            public void handle(ActionEvent event) {
                removeTileListeners();
                listview_polytiles.getItems().clear();
                canvas.resetBlank(canvas.getWidth(), canvas.getHeight());
                updateable.set(true);
            }
        });

        menu_clear_glues.setOnAction(new EventHandler<ActionEvent>() {
            @Override
            public void handle(ActionEvent event) {
                for (Iterator<Glue> iterator = table_gluetable.getItems().iterator(); iterator.hasNext();) {
                    Glue g = iterator.next();
                    iterator.remove();
                    Pair<String, String> p = new Pair<String, String>(g.getGlueA(), g.getGlueB());
                    tc.getGlueFunction().remove(p);
                }
                updateable.set(true);
            }
        });

        menu_export.setOnAction(new EventHandler<ActionEvent>() {
            @Override
            public void handle(ActionEvent event) {
                FileChooser fileChooser = new FileChooser();
                fileChooser.setTitle("Export Tile Configuration...");
                fileChooser.setSelectedExtensionFilter(new FileChooser.ExtensionFilter("XML", "*.xml"));

                File selectedFile = fileChooser.showSaveDialog(null);
                if (selectedFile != null) {
                    try{
                        JAXBContext jaxbContext = JAXBContext.newInstance(TileConfiguration.class);
                        Marshaller marshaller = jaxbContext.createMarshaller();
                        marshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, true);
                        marshaller.marshal(tc, selectedFile);
                    }catch(JAXBException jxb){

                    }
                }

            }
        });

        menu_close.setOnAction(new EventHandler<ActionEvent>() {
            @Override
            public void handle(ActionEvent event) {
                Stage stage  = (Stage) ptAnchorPane.getScene().getWindow();
                // do what you have to do
                stage.close();
            }
        });

        btn_update_assembly.setOnAction(new EventHandler<ActionEvent>() {
            @Override
            public void handle(ActionEvent event) {
                TileSystem ts = current.assembly.getTileSystem();
                ts.getTileTypes().clear();
                ts.getGlueFunction().clear();
                ts.loadTileConfiguration(tc);
                current.removeFrontierFromGrid();
                current.assembly.changeTileSystem(ts);
                updateable.set(false);
            }
        });

    }
}
