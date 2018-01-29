define("galaxy.pages",["exports","mvc/grid/grid-view","libs/jquery/jquery.form","libs/jquery/jstorage","libs/jquery/jquery.wymeditor","libs/jquery/jquery.autocomplete"],function(e,t){"use strict";function a(e){var t,a,i,o;switch(e){case l.ITEM_HISTORY:t="History",a="Histories",i="history",o="History";break;case l.ITEM_DATASET:t="Dataset",a="Datasets",i="dataset",o="HistoryDatasetAssociation";break;case l.ITEM_WORKFLOW:t="Workflow",a="Workflows",i="workflow",o="StoredWorkflow";break;case l.ITEM_PAGE:t="Page",a="Pages",i="page",o="Page";break;case l.ITEM_VISUALIZATION:t="Visualization",a="Visualizations",i="visualization",o="Visualization"}var r="list_"+a.toLowerCase()+"_for_selection";return{singular:t,plural:a,controller:i,iclass:o,list_ajax_url:list_objects_url.replace("LIST_ACTION",r)}}function i(e,t,a){var i=set_accessible_url.replace("ITEM_CONTROLLER",e);$.ajax({type:"POST",url:i,data:{id:t,accessible:"True"},error:function(){alert("Making "+a+" accessible failed")}})}Object.defineProperty(e,"__esModule",{value:!0}),e.default=function(){$(document).ajaxError(function(e,t){var a=t.responseText||t.statusText||"Could not connect to server";return show_modal("Server error",a,{"Ignore error":hide_modal}),!1}),$("[name=page_content]").wymeditor({skin:"galaxy",basePath:editor_base_path,iframeBasePath:iframe_base_path,boxHtml:"<table class='wym_box' width='100%' height='100%'><tr><td><div class='wym_area_top'>"+WYMeditor.TOOLS+"</div></td></tr><tr height='100%'><td><div class='wym_area_main' style='height: 100%;'>"+WYMeditor.IFRAME+WYMeditor.STATUS+"</div></div></td></tr></table>",toolsItems:[{name:"Bold",title:"Strong",css:"wym_tools_strong"},{name:"Italic",title:"Emphasis",css:"wym_tools_emphasis"},{name:"Superscript",title:"Superscript",css:"wym_tools_superscript"},{name:"Subscript",title:"Subscript",css:"wym_tools_subscript"},{name:"InsertOrderedList",title:"Ordered_List",css:"wym_tools_ordered_list"},{name:"InsertUnorderedList",title:"Unordered_List",css:"wym_tools_unordered_list"},{name:"Indent",title:"Indent",css:"wym_tools_indent"},{name:"Outdent",title:"Outdent",css:"wym_tools_outdent"},{name:"Undo",title:"Undo",css:"wym_tools_undo"},{name:"Redo",title:"Redo",css:"wym_tools_redo"},{name:"CreateLink",title:"Link",css:"wym_tools_link"},{name:"Unlink",title:"Unlink",css:"wym_tools_unlink"},{name:"InsertImage",title:"Image",css:"wym_tools_image"},{name:"InsertTable",title:"Table",css:"wym_tools_table"}]});var e=$.wymeditors(0),t=function(t){show_modal("Saving page","progress"),$.ajax({url:save_url,type:"POST",data:{id:page_id,content:e.xhtml(),annotations:JSON.stringify(new Object),_:"true"},success:function(){t()}})};$("#save-button").click(function(){t(function(){hide_modal()})}),$("#close-button").click(function(){window.document.location=page_list_url});var a=$("<div class='galaxy-page-editor-button'><a id='insert-galaxy-link' class='action-button popup' href='#'>Paragraph type</a></div>");$(".wym_area_top").append(a);var i={};$.each(e._options.containersItems,function(t,a){var o=a.name;i[a.title.replace("_"," ")]=function(){e.container(o)}}),make_popupmenu(a,i);var o=$("<div><a id='insert-galaxy-link' class='action-button popup' href='#'>Insert Link to Galaxy Object</a></div>").addClass("galaxy-page-editor-button");$(".wym_area_top").append(o),make_popupmenu(o,{"Insert History Link":function(){e.dialog(l.DIALOG_HISTORY_LINK)},"Insert Dataset Link":function(){e.dialog(l.DIALOG_DATASET_LINK)},"Insert Workflow Link":function(){e.dialog(l.DIALOG_WORKFLOW_LINK)},"Insert Page Link":function(){e.dialog(l.DIALOG_PAGE_LINK)},"Insert Visualization Link":function(){e.dialog(l.DIALOG_VISUALIZATION_LINK)}});var r=$("<div><a id='embed-galaxy-object' class='action-button popup' href='#'>Embed Galaxy Object</a></div>").addClass("galaxy-page-editor-button");$(".wym_area_top").append(r),make_popupmenu(r,{"Embed History":function(){e.dialog(l.DIALOG_EMBED_HISTORY)},"Embed Dataset":function(){e.dialog(l.DIALOG_EMBED_DATASET)},"Embed Workflow":function(){e.dialog(l.DIALOG_EMBED_WORKFLOW)},"Embed Visualization":function(){e.dialog(l.DIALOG_EMBED_VISUALIZATION)}})};var o=function(e){return e&&e.__esModule?e:{default:e}}(t),l={ITEM_HISTORY:"item_history",ITEM_DATASET:"item_dataset",ITEM_WORKFLOW:"item_workflow",ITEM_PAGE:"item_page",ITEM_VISUALIZATION:"item_visualization",DIALOG_HISTORY_LINK:"link_history",DIALOG_DATASET_LINK:"link_dataset",DIALOG_WORKFLOW_LINK:"link_workflow",DIALOG_PAGE_LINK:"link_page",DIALOG_VISUALIZATION_LINK:"link_visualization",DIALOG_EMBED_HISTORY:"embed_history",DIALOG_EMBED_DATASET:"embed_dataset",DIALOG_EMBED_WORKFLOW:"embed_workflow",DIALOG_EMBED_PAGE:"embed_page",DIALOG_EMBED_VISUALIZATION:"embed_visualization"};WYMeditor.editor.prototype.dialog=function(e,t,r){var s=this,n=s.uniqueStamp(),_=s.selected();if(e==WYMeditor.DIALOG_LINK){_&&($(s._options.hrefSelector).val($(_).attr(WYMeditor.HREF)),$(s._options.srcSelector).val($(_).attr(WYMeditor.SRC)),$(s._options.titleSelector).val($(_).attr(WYMeditor.TITLE)),$(s._options.altSelector).val($(_).attr(WYMeditor.ALT)));var d,c;_&&(void 0==(d=$(_).attr("href"))&&(d=""),void 0==(c=$(_).attr("title"))&&(c="")),show_modal("Create Link","<div><div><label id='link_attribute_label'>URL <span style='float: right; font-size: 90%'><a href='#' id='set_link_id'>Create in-page anchor</a></span></label><br><input type='text' class='wym_href' value='"+d+"' size='40' /></div><div><label>Title</label><br><input type='text' class='wym_title' value='"+c+"' size='40' /></div><div>",{"Make link":function(){var e=$(s._options.hrefSelector).val()||"",t=$(".wym_id").val()||"",a=$(s._options.titleSelector).val()||"";if(e||t){s._exec(WYMeditor.CREATE_LINK,n);var i=$("a[href="+n+"]",s._doc.body);i.attr(WYMeditor.HREF,e).attr(WYMeditor.TITLE,a).attr("id",t),0===i.text().indexOf("wym-")&&i.text(a)}hide_modal()},Cancel:function(){hide_modal()}},{},function(){$("#set_link_id").click(function(){$("#link_attribute_label").text("ID/Name");var e=$(".wym_href");e.addClass("wym_id").removeClass("wym_href"),_&&e.val($(_).attr("id")),$(this).remove()})})}if(e==WYMeditor.DIALOG_IMAGE)return s._selected_image&&($(s._options.dialogImageSelector+" "+s._options.srcSelector).val($(s._selected_image).attr(WYMeditor.SRC)),$(s._options.dialogImageSelector+" "+s._options.titleSelector).val($(s._selected_image).attr(WYMeditor.TITLE)),$(s._options.dialogImageSelector+" "+s._options.altSelector).val($(s._selected_image).attr(WYMeditor.ALT))),void show_modal("Image","<div class='row'><label>URL</label><br><input type='text' class='wym_src' value='' size='40' /></div><div class='row'><label>Alt text</label><br><input type='text' class='wym_alt' value='' size='40' /></div><div class='row'><label>Title</label><br><input type='text' class='wym_title' value='' size='40' /></div>",{Insert:function(){var e=$(s._options.srcSelector).val();e.length>0&&(s._exec(WYMeditor.INSERT_IMAGE,n),$("img[src$="+n+"]",s._doc.body).attr(WYMeditor.SRC,e).attr(WYMeditor.TITLE,$(s._options.titleSelector).val()).attr(WYMeditor.ALT,$(s._options.altSelector).val())),hide_modal()},Cancel:function(){hide_modal()}});if(e==WYMeditor.DIALOG_TABLE&&show_modal("Table","<div class='row'><label>Caption</label><br><input type='text' class='wym_caption' value='' size='40' /></div><div class='row'><label>Summary</label><br><input type='text' class='wym_summary' value='' size='40' /></div><div class='row'><label>Number Of Rows<br></label><input type='text' class='wym_rows' value='3' size='3' /></div><div class='row'><label>Number Of Cols<br></label><input type='text' class='wym_cols' value='2' size='3' /></div>",{Insert:function(){var e=$(s._options.rowsSelector).val(),t=$(s._options.colsSelector).val();if(e>0&&t>0){var a=s._doc.createElement(WYMeditor.TABLE),i=null,o=$(s._options.captionSelector).val();a.createCaption().innerHTML=o;for(var l=0;l<e;l++){i=a.insertRow(l);for(var r=0;r<t;r++)i.insertCell(r)}$(a).attr("summary",$(s._options.summarySelector).val());var n=$(s.findUp(s.container(),WYMeditor.MAIN_CONTAINERS)).get(0);n&&n.parentNode?$(n).after(a):$(s._doc.body).append(a)}hide_modal()},Cancel:function(){hide_modal()}}),e==l.DIALOG_HISTORY_LINK||e==l.DIALOG_DATASET_LINK||e==l.DIALOG_WORKFLOW_LINK||e==l.DIALOG_PAGE_LINK||e==l.DIALOG_VISUALIZATION_LINK){switch(e){case l.DIALOG_HISTORY_LINK:I=a(l.ITEM_HISTORY);break;case l.DIALOG_DATASET_LINK:I=a(l.ITEM_DATASET);break;case l.DIALOG_WORKFLOW_LINK:I=a(l.ITEM_WORKFLOW);break;case l.DIALOG_PAGE_LINK:I=a(l.ITEM_PAGE);break;case l.DIALOG_VISUALIZATION_LINK:I=a(l.ITEM_VISUALIZATION)}u=new o.default({url_base:I.list_ajax_url,dict_format:!0,embedded:!0});Galaxy.modal.show({title:"Insert Link to "+I.singular,body:$("<div/>").append(u.$el).append($("<div/>").append('<input id="make-importable" type="checkbox" checked/>').append("Make the selected "+I.plural.toLowerCase()+" accessible so that they can viewed by everyone.")),closing_events:!0,buttons:{Insert:function(){var e=!1;null!=$("#make-importable:checked").val()&&(e=!0);new Array;u.$("input[name=id]:checked").each(function(){var t=$(this).val();e&&i(I.controller,t,I.singular);var a=(get_name_and_link_url+t).replace("ITEM_CONTROLLER",I.controller);$.getJSON(a,function(e){s._exec(WYMeditor.CREATE_LINK,n);var a=$("a[href="+n+"]",s._doc.body).text();""==a||a==n?s.insert("<a href='"+e.link+"'>"+I.singular+" '"+e.name+"'</a>"):$("a[href="+n+"]",s._doc.body).attr(WYMeditor.HREF,e.link).attr(WYMeditor.TITLE,I.singular+t)})}),Galaxy.modal.hide()},Close:function(){Galaxy.modal.hide()}}})}if(e==l.DIALOG_EMBED_HISTORY||e==l.DIALOG_EMBED_DATASET||e==l.DIALOG_EMBED_WORKFLOW||e==l.DIALOG_EMBED_PAGE||e==l.DIALOG_EMBED_VISUALIZATION){var I;switch(e){case l.DIALOG_EMBED_HISTORY:I=a(l.ITEM_HISTORY);break;case l.DIALOG_EMBED_DATASET:I=a(l.ITEM_DATASET);break;case l.DIALOG_EMBED_WORKFLOW:I=a(l.ITEM_WORKFLOW);break;case l.DIALOG_EMBED_PAGE:I=a(l.ITEM_PAGE);break;case l.DIALOG_EMBED_VISUALIZATION:I=a(l.ITEM_VISUALIZATION)}var u=new o.default({url_base:I.list_ajax_url,dict_format:!0,embedded:!0});Galaxy.modal.show({title:"Insert Link to "+I.singular,body:$("<div/>").append(u.$el).append($("<div/>").append('<input id="make-importable" type="checkbox" checked/>').append("Make the selected "+I.plural.toLowerCase()+" accessible so that they can viewed by everyone.")),closing_events:!0,buttons:{Embed:function(){var e=!1;null!=$("#make-importable:checked").val()&&(e=!0),u.$("input[name=id]:checked").each(function(){var t=$(this).val(),a=$("label[for='"+t+"']:first").text();e&&i(I.controller,t,I.singular);var o=["<div id='",I.iclass+"-"+t,"' class='embedded-item ",I.singular.toLowerCase()," placeholder'>","<p class='title'>","Embedded Galaxy ",I.singular," '",a,"'","</p>","<p class='content'>","[Do not edit this block; Galaxy will fill it in with the annotated ",I.singular.toLowerCase()," when it is displayed.]","</p>","</div>"].join("");s.insert(o)}),Galaxy.modal.hide()},Close:function(){Galaxy.modal.hide()}}})}}});
//# sourceMappingURL=../maps/galaxy.pages.js.map