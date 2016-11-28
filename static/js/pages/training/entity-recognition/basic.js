YPet = new Backbone.Marionette.Application({
    region: "#basics"
});

var load_initial_data;

(load_initial_data = function() {
    var p = new Paragraph({'text': $('#basics').text().trim()});
    YPet.AnnotationTypes = new AnnotationTypeList([{name: 'Disease', color: '#d1f3ff'}]);
    // YPet.addRegions({'basics': '#basics'});
    YPet.getRegion().show(new WordCollectionView({collection: p.get('words')}))
})();

YPet.start();


$('a').on('click', function(evt) {
    if ($(this).hasClass('disabled')) {
        evt.preventDefault();
    }
});

$('#next-container').popover({
    trigger: "hover",
    title: "Keep Going!",
    content: "Please complete the task before moving forward.",
    placement: "top",
    container: "body"
});

var step_function = function() {
    switch (step_idx) {
        case 0:
        case 1:
            $biomedical_el.popover({
                container: 'body',
                html: true,
                animation: false,
                content: function() {
                    if (step_idx == 0) {
                        return 'Highlight words by clicking on them. Highlight the word \“<span class="user_annotated">Biomedical</span>\”.';
                    } else {
                        return 'Remove the highlight by clicking again. Un-highlight the word \“<span class="user_annotated">Biomedical</span>\”.';
                    }
                },
                placement: 'bottom'
            });
            $biomedical_el.popover('show');
            break;

        case 2:
            $biomedical_el.popover('hide');
            $because_el.popover({
                container: 'body',
                html: true,
                animation: false,
                content: 'Highlight blocks of text by holding down your mouse button until you reach the end of the text you wish to highlight. Highlight the phrase \“<span class="user_annotated">because I am ready to help</span>\”.',
                placement: 'top'
            });
            $because_el.popover('show');
            break;

        case 3:
            $biomedical_el.popover('hide');
            $because_el.popover('hide');
            $next_el.popover({
                container: 'body',
                html: true,
                animation: false,
                content: 'Click Next to continue!',
                placement: 'top'
            });
            $next_el.popover('show');
            $next_el.attr('disabled', false).removeClass('disabled');
            $('#next-container').popover('disable');
            break;
    }
    ;

};

var step_idx = 0;
$('body').ready(function() {
    $biomedical_el = $('p.paragraph span:nth(0)');
    $because_el = $('p.paragraph  span:nth(56)');
    $help_el = $('p.paragraph span:nth(61)');
    same_line = $because_el.position().top == $help_el.position().top;
    $next_el = $('#next');

    step_function();
});

YPet.getRegion().currentView.collection.parentDocument.get('annotations').on('add', function(model, collection) {
    var model_json = model.toJSON();

    /* Unselect wrong answers */
    if (!_.contains(['Biomedical', 'because I am ready to help'], model_json.text) && model.get('words').first()) {

        setTimeout(function() {
            var first_word = model.get('words').first();
            first_word.trigger('unclick');
        }, 50);

        collection.drawAnnotations(model);
    }

    if (model_json.text == "Biomedical" && collection.length == 1) {
        step_idx = 1;
        step_function();
    }
    if (model_json.text == "because I am ready to help" && collection.length == 1) {
        step_idx = 3;
        step_function();
    }
});

YPet.getRegion().currentView.collection.parentDocument.get('annotations').on('remove', function(model, collection) {
    var model_json = model.toJSON();
    if (model_json.text == "Biomedical" && collection.length == 0) {
        step_idx = 2;
        step_function();
    }
});
