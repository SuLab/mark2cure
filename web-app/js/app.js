define(['marionette', 'vent',
        //-- Models
        'models/Document', 'models/User',

        //-- Collections
        'collections/DocumentList',

        //-- Views
        'views/Header/Index', 'views/Footer/Index',
        'views/Main/Game/Index', 'views/Main/Library/Index',
        'views/Analytics/Network',
        'views/modals/FirstRun', 'views/modals/Complete', 'views/modals/Message', 'views/modals/Survey',
        'views/modals/Settings', 'views/modals/Instructions',

        //-- Utils
        'extensions/ModalRegion', 'extensions/FadeRegion', 'extensions/SlideRegion'],

   function(marionette, vent,
            //-- Models
            Document, User,

            //-- Collections
            DocumentList,

            //-- Views
            Header, Footer,
            Game, Library,
            Network,
            FirstRun, Complete, Message, Survey,
            Settings, Instructions,

            //--Utils
            ModalRegion, FadeRegion, SlideRegion) {
    'use strict';

    var app = new marionette.Application(),
        viewOptions = {
          collection : new DocumentList(),
          user : User
        };

    app.addRegions({
      header    : '#header',
      analytics : SlideRegion.extend({el: "#analytics"}),
      main      : FadeRegion.extend({el: "#content"}),
      footer    : '#footer',
      modal     : ModalRegion,
    });

    app.addInitializer(function() {
      viewOptions.collection.fetch({success : function() {
        app.header.show(  new Header(viewOptions)  );
        app.main.show(    new Library(viewOptions) );
        app.footer.show(  new Footer(viewOptions)  );

        if( viewOptions.user.get('first_run') ) {
          app.main.close();
          app.footer.close();

          //-- This set the modal so it can't be closed
          viewOptions.static = true;
          app.modal.show(   new FirstRun(viewOptions) )
        }
      }});

    });

    vent.on('navigate', function(param) {
      if( viewOptions.user.get('first_run') ) {
        // app.header.close();
      }

      var found_model = viewOptions.collection.findWhere({id: Number(param) });
      if(found_model) {
        viewOptions['model'] = found_model;
        app.main.show( new Game(viewOptions) );
      } else {
        Backbone.history.navigate( '#/library' );
        vent.trigger('library', {});
      }
    });

    vent.on('userchanged', function() {
      viewOptions.collection.fetch();
    });

    vent.on('navigate:library', function(obj) {
      var quest = obj['quest']

      if(quest && quest.trim().length) {
        $.ajax({
          async: false,
          url: '/api/v1/quest/'+quest.trim(),
          success: function(jsonData) {
            var docs = []
            _.each(jsonData.objects, function(v) {
              var doc = new Document({id: v.document_id});
              doc.fetch();
              docs.push(doc);
            });
            viewOptions.collection = new DocumentList(docs);
            app.main.show( new Library(viewOptions) );

            if(docs.length==0) {
              Backbone.history.navigate( '#/library' );
              vent.trigger('library', {});
            }

          },

          error: function() {
            Backbone.history.navigate( '#/library' );
            vent.trigger('library', {});
          }
        });

      } else {
        app.main.show( new Library(viewOptions) );
      }
    });

    vent.on('navigate:analytics', function(obj) {
      if(!obj['toggle'] || $('#network').length ) {
        app.analytics.close();
      } else {
        app.analytics.show( new Network(obj) );
      }
    });

    //
    //-- Display Modals
    //
    vent.on('navigate:message', function() {
      app.modal.show( new Message(viewOptions) );
    });

    vent.on('navigate:survey', function() {
      app.modal.show( new Survey(viewOptions) );
    });

    vent.on('navigate:settings', function() {
      app.modal.show(   new Settings(viewOptions) )
    });

    vent.on('navigate:instructions', function() {
      app.modal.show(   new Instructions(viewOptions) )
    });

    //-- Oneoff modals really...
    vent.on('modal:show_complete', function() {
      viewOptions['swap'] = true;
      app.modal.show( new Complete(viewOptions) );
    });

    vent.on('modal:welcome', function() {
      app.modal.show( new FirstRun(viewOptions) );
    });

    vent.on('modal:close', function() {
      Backbone.history.navigate('#');
      $('#modal').modal('hide');
    });

    return app;

  }
);
