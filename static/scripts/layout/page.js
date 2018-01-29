define("layout/page",["exports","layout/masthead","layout/panel","mvc/ui/ui-modal","utils/utils"],function(e,i,t,s,a){"use strict";function n(e){return e&&e.__esModule?e:{default:e}}Object.defineProperty(e,"__esModule",{value:!0});var o=n(i),l=n(t),r=n(s),d=n(a),c=Backbone.View.extend({el:"body",className:"full-content",_panelids:["left","right"],initialize:function(e){var i=this;this.config=_.defaults(e.config||{},{message_box_visible:!1,message_box_content:"",message_box_class:"info",show_inactivity_warning:!1,inactivity_box_content:""}),Galaxy.modal=this.modal=new r.default.View,Galaxy.display=this.display=function(e){e.title?(d.default.setWindowTitle(e.title),e.allow_title_display=!1):(d.default.setWindowTitle(),e.allow_title_display=!0),i.center.display(e)},Galaxy.router=this.router=e.Router&&new e.Router(i,e),this.masthead=new o.default.View(this.config),this.center=new l.default.CenterPanel,this.$el.attr("scroll","no"),this.$el.html(this._template()),this.$("#masthead").replaceWith(this.masthead.$el),this.$("#center").append(this.center.$el),this.$el.append(this.masthead.frame.$el),this.$el.append(this.modal.$el),this.$messagebox=this.$("#messagebox"),this.$inactivebox=this.$("#inactivebox"),this.panels={},_.each(this._panelids,function(t){var s=t.charAt(0).toUpperCase()+t.slice(1),a=e[s];if(a){var n=new a(i,e);i[n.toString()]=n,i.panels[t]=new l.default.SidePanel({id:t,el:i.$("#"+t),view:n})}}),this.render(),this.router&&Backbone.history.start({root:Galaxy.root,pushState:!0})},render:function(){return $(".select2-hidden-accessible").remove(),this.masthead.render(),this.renderMessageBox(),this.renderInactivityBox(),this.renderPanels(),this._checkCommunicationServerOnline(),this},renderMessageBox:function(){if(this.config.message_box_visible){var e=this.config.message_box_content||"",i=this.config.message_box_class||"info";this.$el.addClass("has-message-box"),this.$messagebox.attr("class","panel-"+i+"-message").html(e).toggle(!!e).show()}else this.$el.removeClass("has-message-box"),this.$messagebox.hide();return this},renderInactivityBox:function(){if(this.config.show_inactivity_warning){var e=this.config.inactivity_box_content||"",i=$("<a/>").attr("href",Galaxy.root+"user/resend_verification").text("Resend verification");this.$el.addClass("has-inactivity-box"),this.$inactivebox.html(e+" ").append(i).toggle(!!e).show()}else this.$el.removeClass("has-inactivity-box"),this.$inactivebox.hide();return this},renderPanels:function(){var e=this;return _.each(this._panelids,function(i){var t=e.panels[i];t?t.render():(e.$("#center").css(i,0),e.$("#"+i).hide())}),this},_template:function(){return['<div id="everything">','<div id="background"/>','<div id="masthead"/>','<div id="messagebox"/>','<div id="inactivebox" class="panel-warning-message" />','<div id="left" />','<div id="center" />','<div id="right" />',"</div>",'<div id="dd-helper" />'].join("")},toString:function(){return"PageLayoutView"},_checkCommunicationServerOnline:function(){var e=window.Galaxy.config.communication_server_host,i=window.Galaxy.config.communication_server_port,t=window.Galaxy.user.attributes.preferences,s=$("#show-chat-online");t&&-1!=["1","true"].indexOf(t.communication_server)?$.ajax({url:e+":"+i}).success(function(e){null!==window.Galaxy.user.id&&"hidden"===s.css("visibility")&&s.css("visibility","visible")}).error(function(e){s.css("visibility","hidden")}):s.css("visibility","hidden")}});e.default={View:c}});
//# sourceMappingURL=../../maps/layout/page.js.map